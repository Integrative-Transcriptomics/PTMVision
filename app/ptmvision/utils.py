import pandas as pd
import numpy as np
import scipy as sc
import pyteomics.mass.unimod as unimod
from Bio.PDB import PDBParser
from Bio import SeqIO
from io import StringIO
from pathlib import Path
from psm_utils.io import read_file
from tempfile import NamedTemporaryFile
from datetime import datetime
import re, urllib.request, requests, brotli, base64, os

BASEPATH = "./app/ptmvision" # Use this for local development.
# BASEPATH = "/app/ptmvision" # Use this for deployment.
UNIMOD_MAPPER = unimod.Unimod("sqlite:///unimod.db") # Use this for local development.
# UNIMOD_MAPPER = unimod.Unimod( ) # Use this for deployment. TODO: Seems to be broken.
TOLERANCE = 0.001 # Mass tolerance when matching masses to unimod IDs.
EXCLUDE_CLASSES = [ ] # Exclude PTMs with these Unimod classifications from the results.
PDBPARSER = PDBParser(PERMISSIVE=False)
DEBUG = os.getenv("DEBUG")
"""
TODO: Refactor into reader classes, one for each file format.
"""

def get_distance_matrix(structure):
    residue_count = 0
    for residue in enumerate(structure.get_residues()):
        residue_count += 1
    coords = np.zeros(shape=(residue_count, 3))

    for i, residue in enumerate(structure.get_residues()):
        if residue.get_resname() == "GLY":
            _ = "CA"
        else:
            _ = "CB"
        for atom in residue:
            if atom.get_name() == _:
                coords[i] = atom.get_coord()
                break
    residue_distances = sc.spatial.distance_matrix(coords, coords)
    return residue_distances


def get_contacts(dist_matrix, threshold):
    contacts = {}
    for pairwise_contact in np.argwhere(dist_matrix < threshold):
        if pairwise_contact[0] != pairwise_contact[1]:
            contacts.setdefault(
                int(pairwise_contact[0]), []
            ).append( int( pairwise_contact[1] ) )
    return contacts


def get_structure(uniprot_id):
    url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb".format(uniprot_id)
    structure = None
    try:
        response = urllib.request.urlopen(url)
        data = response.read()
        pdb_text = data.decode("utf-8")
        structure = PDBPARSER.get_structure(uniprot_id, StringIO(pdb_text))
    except:
        raise Exception("Alphafold structure not available.")
    return structure, pdb_text


def parse_structure(structure_string):
    structure = None
    # with warnings.filterwarnings("error"):
    try:
        structure = PDBPARSER.get_structure("custom_pdb", StringIO(structure_string))
    except Exception as e:
        raise Exception("Failed to parse provided .pdb structure: " + str(e))
    finally:
        return structure, structure_string


def query_uniprot_for_fasta(uniprot_ids):
    """
    request sequences from uniprot via REST interface from uniprot IDs
    """
    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28"
    url_separator = "%20OR%20"
    id_list = [
        "%28accession%3A{}".format(uniprot_id + "%29") for uniprot_id in uniprot_ids
    ]

    fasta_text = ""

    # split into batches to avoid bad gateway
    id_lists = [id_list[x : x + 100] for x in range(0, len(id_list), 100)]
    for id_chunk in id_lists:
        id_chunk = url_separator.join(id_chunk) + "%29"
        query = url + id_chunk
        response = requests.get(query)

        if response.status_code == 200:
            response_text = "".join(response.text)
            fasta_text += response_text

    fasta = SeqIO.parse(StringIO(fasta_text), "fasta")
    fasta_dict = {}
    for rec in fasta:
        fasta_dict[rec.id.split("|")[1]] = str(rec.seq)
    return fasta_dict


def get_protein_sequence(fasta_dict, uniprot_id):
    try:
        return fasta_dict[uniprot_id]
    except KeyError:
        return ""


def get_sequence_from_structure(structure):
    sequence = []
    for residue in structure.get_residues():
        sequence.append(residue.resname)
    return sequence


def get_peptide_position_in_protein(peptide, protein):
    if protein == "":
        return -1
    peptide_sequence = (
        re.sub(r"\[(.*?)\]", "", peptide)
        .replace("[", "")
        .replace("]", "")
        .replace("-", "")
    )
    return int(protein.find(peptide_sequence))


def get_amino_acid_from_protein(position, uniprot_id, fasta_dict):
    for record in fasta_dict.keys():
        if uniprot_id in record:
            return fasta_dict[record][position - 1]  # positions are 1 indexed
    return "N/A"


def gather_mass_shift_msfragger(row):
    """
    gather mass shifts from MSFragger PSM row
    1) try to parse from modifications tuple
    2) then from Observed Modifications	(unannotated mass shifts)
    3) then from unimod id mapping
    """
    if not pd.isna(row["modification"]):
        mod = row["modification"][1]
        try:
            mass_shift = float(mod)
            return mass_shift
        except ValueError:
            pass

    if not pd.isna(row["Observed Modifications"]):
        if "Unannotated mass-shift" in row["Observed Modifications"]:
            match = match = re.search(
                r"Unannotated mass-shift (-?\d+(\.\d+)?)", row["Observed Modifications"]
            )
            if match:
                return float(match.group(1))
            else:
                return None

    if not pd.isna(row["modification_unimod_id"]):
        # just check one representative mass shift
        mass_shift = map_id_to_mass(str(row["modification_unimod_id"].split(" or ")[0]))
        return mass_shift

    return


def map_modification_string_to_unimod_id(modification_string):
    # gets either name ("Oxidation") or string ("+57.") and returns unimod id(s)
    if modification_string[0] in ["+", "-"]:
        return map_mass_to_unimod_id(float(modification_string))
    else:
        # check if valid unimod name
        query_result = UNIMOD_MAPPER.by_name(modification_string, strict=False)
        if query_result:
            return str(query_result.id)
        else:
            return ""


def map_mass_to_unimod_id(mass):
    global TOLERANCE
    query_results = UNIMOD_MAPPER.session.query(unimod.Modification).filter(unimod.Modification.monoisotopic_mass.between(mass-TOLERANCE, mass+TOLERANCE)).all()
    mod_ids = [str(result.id) for result in query_results]
    if len(mod_ids) == 0:
        return ""
    elif len(mod_ids) > 1:
        return " or ".join(mod_ids)
    return str(mod_ids[0])


def map_mass_to_unimod_names(mass):
    global TOLERANCE
    query_results = UNIMOD_MAPPER.session.query(unimod.Modification).filter(unimod.Modification.monoisotopic_mass.between(mass-TOLERANCE, mass+TOLERANCE)).all()
    mod_names = [str(result.code_name) for result in query_results]
    if len(mod_names) == 0:
        return ""
    elif len(mod_names) > 1:
        return " or ".join(mod_names)
    return str(mod_names[0])

def map_unimod_name_to_mass(unimod_name):
    query_results = UNIMOD_MAPPER.session.query(unimod.Modification).filter(unimod.Modification.ex_code_name == unimod_name).all()
    mod_masses = [result.monoisotopic_mass for result in query_results]
    if len(mod_masses) == 0:
        return None
    elif len(mod_masses) > 1:
        return " or ".join(mod_masses)
    return mod_masses[0]


def map_ids_to_unimod_name(ids):
    names = []
    for id in ids.split(" or "):
        name = map_id_to_unimod_name(id)
        if name:
            names.append(name)
    if len(names) > 0:
        return " or ".join(names)
    return None


def map_id_to_unimod_name(id):
    query_results = UNIMOD_MAPPER.session.query(unimod.Modification).filter(unimod.Modification.id == id).all()
    mod_names = [result.ex_code_name if result.ex_code_name != "" else result.code_name for result in query_results  ]
    if len(mod_names) == 0:
        return None
    elif len(mod_names) > 1:
        return " or ".join(mod_names)
    return mod_names[0]


def map_id_to_mass(id):
    query_results = UNIMOD_MAPPER.session.query(unimod.Modification).filter(unimod.Modification.id == id).all()
    mod_masses = [result.monoisotopic_mass for result in query_results]
    if len(mod_masses) == 0:
        return None
    elif len(mod_masses) > 1:
        return " or ".join(mod_masses)
    return mod_masses[0]


def map_unimod_description_to_unimod_id(mod):
    query_results = UNIMOD_MAPPER.session.query(unimod.Modification).filter(unimod.Modification.full_name == mod).all()
    mod_ids = [str(result.id) for result in query_results]
    if len(mod_ids) == 0:
        return None
    elif len(mod_ids) > 1:
        return " or ".join(mod_ids)
    return mod_ids[0]


def mass_shift_from_id(unimod_id, unimod_db):
    """
    Get mass shift from unimod ID
    """
    try:
        mod = unimod_db.get(unimod_id)
    except KeyError:  # deprecated unimod ID :(
        return None
    return mod.monoisotopic_mass


def unimod_ids_to_unimod_names(x):
    """
    map list of unimod ids to list of unimod names
    """
    if x is None:
        return None
    names = []
    for x in x.split(" or "):
        names.append(map_id_to_unimod_name(x))
    names = [x for x in names if x is not None]
    if len(names) > 0:
        return " or ".join(names)
    return


def classification_from_id(unimod_id, amino_acid, unimod_db):
    """
    Get amino acid specific classification of modification from unimod ID
    """
    if unimod_id == "": # If the mass shift is not annotated in Unimod.
        return "Unannotated mass shift"
    elif type(unimod_id) == str and " or " in unimod_id: # If the mass shift is ambiguous, i.e. matches multiple Unimod IDs.
        return "Multiple"
    else: # If the mass shift is annotated in Unimod, try to match.
        try:
            mod = unimod_db.get(int(unimod_id))
        except KeyError:  # In case of deprecated Unimod ID :(
            return "Deprecated unimod entry"
    for specification in mod.specificities:
        if specification.amino_acid == amino_acid:
            classification = specification.classification.classification
            return classification
    return "Non-standard residue" # If the mass shift is not annotated for the amino acid.

def get_classification(unimod_id, residue, unimod_db):
    if unimod_id is None:
        return "Unannotated mass shift"
    elif type(unimod_id) == str and " or " in unimod_id:
        return "Multiple"
    else:
        classification = classification_from_id(unimod_id, residue, unimod_db)
        return classification


def extract_mods_from_proforma(peptidoform):
    """
    ENYC[+57.0215]NNVMM[+15.994915]K/2	-> [(3, +57.0215), (8, +15.993915)] (zero indexed!), try to map to unimod names
    AADM[Oxidation]TGADIEAMTR/2 -> [(3, Oxidation)]
    """
    mods = []
    while peptidoform.count("[") > 0:
        # get first match
        match = re.search(r"\[.+?\]", peptidoform)
        if match:
            # get modified residue index, append index and mass shift to list
            modified_residue_index = int(match.start()) - 1
            mod = match.group(0).replace("[", "").replace("]", "")
            peptidoform = peptidoform[: match.start()] + peptidoform[match.end() :]
            mods.append((modified_residue_index, mod))
        else:
            return None
    return mods


def parse_protein_string(x):
    """
    Try to parse uniprot accession ID from protein string
    """
    uniprot_accession_pattern = re.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )
    accession_candidates = [x] + x.split("|")
    for candidate in accession_candidates:
        match = uniprot_accession_pattern.match(candidate)
        if match:
            return match.group(0)  # probably accession number
    return None


def msfragger_to_unimod_id(x):
    """
    map mass shift or description to unimod id
    """
    # (7, 15.9949)
    # (3, Unannotated mass-shift -116.0590)

    # Try mass shift
    mod = x[1]
    try:
        mass_shift = float(mod)
        return map_mass_to_unimod_id(mass_shift)

    # then try description
    except ValueError:
        if "Unannotated mass-shift" in mod:
            mass_shift = float(
                mod.split("Unannotated mass-shift ")[1].split(", Mod2:")[0]
            )
            # try to map to unimod again (will probably not work since PTMShepherd didnt find a match either)
            return map_mass_to_unimod_id(mass_shift)
        unimod_ids = []
        for description in mod.split("/"):
            unimod_ids.append(map_unimod_description_to_unimod_id(description))
        unimod_ids = [x for x in unimod_ids if x is not None]
        if len(unimod_ids) == 0:
            return None
        return " or ".join(unimod_ids)


def get_display_name(row):
    # if unimod name is available and unambiguous, use it as display name
    if row["modification_unimod_name"] and " or " not in row["modification_unimod_name"]:
            return row["modification_unimod_name"]
    # otherwise use the original label
    else:
        return row["modification"][1]


def get_mass_shift(mod):
    # gets either name ("Oxidation") or string ("+57.") and returns mass shift
    if mod[0] in ["+", "-"]:
        mass_shift = float(mod)
    else:
        # try to get mass over unimod mapper
        mass_shift = map_unimod_name_to_mass(mod)
    return mass_shift

def get_candidates_from_mass_shift(mass_shift):
    unimod_names = map_mass_to_unimod_names(mass_shift)
    if " or " in unimod_names:
        return unimod_names
    return

def get_AA(peptide, position):
    return peptide[position]

def PSMList_to_mod_df(psm_list):
    """
    parse list of PSMs and return pandas dataframe with [protein, position, modification_unimod_id, modification_unimod_name]
    # TODO: refactor
    """
    peptidoforms = [psm["peptidoform"].proforma for psm in psm_list]
    proteins = [psm.protein_list[0] for psm in psm_list]
    peptide_sequences = [psm["peptidoform"].sequence for psm in psm_list]

    df = pd.DataFrame(
        zip(peptidoforms, proteins, peptide_sequences),
        columns=["peptidoform", "protein", "peptide"],
    )

    df["modification"] = df["peptidoform"].apply(
        lambda x: extract_mods_from_proforma(x)
    )

    # one mod per line
    df = df.explode("modification")

    df["uniprot_id"] = df["protein"].apply(lambda x: parse_protein_string(x))
    df = df.dropna(subset=["uniprot_id"])

    # query uniprot for protein sequence
    fasta = query_uniprot_for_fasta(df["uniprot_id"].unique().tolist())

    df["protein_sequence"] = df["uniprot_id"].apply(
        lambda x: get_protein_sequence(fasta, x)
    )
    df = df[df["protein_sequence"] != ""]

    df["mod_position_in_peptide"] = df["modification"].apply(lambda x: x[0])

    #get all unique modifications
    unique_mods = df[["modification"]].drop_duplicates()

    #look unimod id up in the db, store them in column
    unique_mods["modification_unimod_id"] = unique_mods["modification"].apply(
        lambda x
        : map_modification_string_to_unimod_id(x[1])
    )

    unique_mods["mass_shift"] = unique_mods["modification"].apply(lambda x: get_mass_shift(x[1]))

    #merge
    df = df.merge(unique_mods, on=["modification"], how="left")

    #get all unique unimod ids
    unique_mods = df[["modification_unimod_id"]].drop_duplicates()

    #look unimod classification up in the db, store them in column
    unique_mods["modification_unimod_name"] = unique_mods["modification_unimod_id"].apply(
        lambda x: map_ids_to_unimod_name(x)
    )

    #merge
    df = df.merge(unique_mods, on=["modification_unimod_id"], how="left")

    df["display_name"] = df.apply(
        lambda x: get_display_name(x), axis=1
    )  # Unimod Name if available, otherwise "unannotated mass shift: (mass shift)"

    df["peptide_position"] = df.apply(
        lambda x: get_peptide_position_in_protein(x["peptide"], x["protein_sequence"]),
        axis=1,
    )

    # get unique unimod - amino acid combinations
    df["AA"] = df.apply(lambda x: get_AA(x["peptide"], x["mod_position_in_peptide"]), axis = 1)

    unique_mod_aas = df[["modification_unimod_id", "AA"]].drop_duplicates()
    unique_mod_aas["classification"] = unique_mod_aas.apply(
        lambda x: classification_from_id(
            x["modification_unimod_id"], x["AA"], UNIMOD_MAPPER
        ),
        axis=1,
    )

    df = df.merge(unique_mod_aas, on=["modification_unimod_id", "AA"], how="left")

    df["position"] = df["mod_position_in_peptide"] + df["peptide_position"]

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    df = df[
        [
            "uniprot_id",
            "position",
            "modification_unimod_name",
            "modification_unimod_id",
            "classification",
            "mass_shift",
            "display_name",
        ]
    ]

    df = df.drop_duplicates()

    return df


def parse_assigned_mods_msfragger(mods):
    "7M(15.9949) => (7, 15.9949), N-Term(42.0106) => (0, Unannotated mass-shift 42.0106)"
    mod_list = []
    for mod in mods.split(","):
        if "N-term" in mod:
            pos = 0
        else:
            pos = int(re.search("\d+", mod).group(0)) - 1
        mass = re.search(r"\(.*?\)", mod)[0]
        mass = float(mass.replace("(", "").replace(")", ""))
        mod_list.append((pos, "Unannotated mass-shift " + str(mass)))

    return mod_list


def parse_observed_mods_msfragger(localisation, mod):
    "Mod1: Methylation (PeakApex: 14.0144, Theoretical: 14.0157) + TVkEAEEAAK => (2, Methylation)"
    match = re.search(r"[a-z]", localisation)
    if match:
        position = match.start()
    else:
        return

    mod_match = re.search(r"Mod1: (.*?)\(", mod)
    if mod_match:
        mod = mod_match.group(1).strip().split(", Mod2")[0]
    else:
        return
    return [(position, mod)]


def parse_msfragger_mods(msfragger_psm_row):
    """
    Parse modifications from MSFragger PSM.

    Variable mods are in column "Assigned Modification" (e.g. "7M(15.9949)" or "N-term(42.0106)").
    Mass shifts are in column "Observed Modification" (e.g. Mod1: Methylation (PeakApex: 14.0144, Theoretical: 14.0157)).
    Their localization is in column "MSFragger Localization" (e.g. TVkEAEEAAK).
    Combine into one list of tuples => [(7, 15.9949), ()].

    :param msfragger_psm_row: A row from the MSFragger PSM DataFrame.
    :return: A list of tuples representing the modifications.
    """
    mods = []

    assigned_mods = msfragger_psm_row.get("Assigned Modifications")
    if pd.notna(assigned_mods):
        var_mods = parse_assigned_mods_msfragger(assigned_mods)
        if var_mods:
            mods.extend(var_mods)

    localization = msfragger_psm_row.get("MSFragger Localization")
    observed_mods = msfragger_psm_row.get("Observed Modifications")
    if pd.notna(localization):
        new_mods = parse_observed_mods_msfragger(localization, observed_mods)
        if new_mods:
            mods.extend(new_mods)

    return mods


def read_msfragger(file):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'msfragger' ... ", end = '', flush = True)
    # read file and pick relevant rows and columns
    pept_mods = pd.read_csv(file, sep="\t")
    pept_mods = pept_mods[
        [
            "Modified Peptide",
            "Peptide",
            "Protein ID",
            "Assigned Modifications",
            "Observed Modifications",
            "MSFragger Localization",
            "Protein Start",
        ]
    ]

    pept_mods = pept_mods.drop_duplicates()

    # filter out unmodified peptidoforms
    pept_mods = pept_mods[
        ~(
            pept_mods["Assigned Modifications"].isna()
            & pept_mods["MSFragger Localization"].isna()
        )
    ]

    # filter out ambigously matched modifications
    pept_mods = pept_mods[pept_mods["Observed Modifications"].fillna("").str.count(";") == 0]

    # filter out combinations of modifications
    pept_mods = pept_mods[pept_mods["Observed Modifications"].fillna("").str.count("Mod2") == 0]

    # filter out ambiguously localized modifications
    pept_mods = pept_mods[
        (
            pept_mods["MSFragger Localization"].apply(
                lambda x: sum(1 for c in str(x) if c.islower())
            )
            == 1 & ~pept_mods["Observed Modifications"].isna()
        )
        | (pept_mods["MSFragger Localization"].isna())
    ]

    # parse modification columns
    pept_mods["modification"] = pept_mods.apply(
        lambda x: parse_msfragger_mods(x), axis=1
    )
    pept_mods = pept_mods.explode("modification")

    # map modification to unimod id
    pept_mods = pept_mods.dropna(subset=["modification"])
    pept_mods["modification_unimod_id"] = pept_mods["modification"].apply(
        lambda x: msfragger_to_unimod_id(x)
    )

    # map modification onto protein position
    # Protein Start is 1 indexed, modification position is 0 indexed, so cancels out the -1
    pept_mods["position"] = pept_mods["Protein Start"] + pept_mods[
        "modification"
    ].apply(lambda x: int(x[0]))

    # parse protein column to accession ID
    pept_mods["uniprot_id"] = pept_mods["Protein ID"].apply(
        lambda x: parse_protein_string(x)
    )
    pept_mods = pept_mods.dropna(subset=["uniprot_id"])

    # map unimod ids to unimod names
    pept_mods["modification_unimod_name"] = pept_mods["modification_unimod_id"].apply(
        lambda x: unimod_ids_to_unimod_names(x)
    )

    # gather mass shift from different columns
    pept_mods["mass_shift"] = pept_mods.apply(
        lambda x: gather_mass_shift_msfragger(x), axis=1
    )

    # get amino acid at modified position and classification of modification
    # TODO: add N-Term Support (currently encoded as pos = 0)
    pept_mods["residue"] = pept_mods.apply(
        lambda x: x["Peptide"][int(x["modification"][0])], axis=1
    )
    pept_mods["classification"] = pept_mods.apply(
        lambda x: get_classification(
            x["modification_unimod_id"], x["residue"], UNIMOD_MAPPER
        ),
        axis=1,
    )

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    pept_mods = pept_mods[~pept_mods["classification"].isin(EXCLUDE_CLASSES)]

    # get display name (unimod name if unambiguous, otherwise "unannotated mass shift: (mass shift)"
    pept_mods["display_name"] = pept_mods.apply(lambda x: get_display_name(x), axis=1)

    # for ambiguous mass shifts, store all candidates
    pept_mods["candidates"] = pept_mods.apply(
        lambda row: row["modification_unimod_name"] if "Multiple" in row["classification"] else None,
        axis=1
    )
    return pept_mods


def read_ionbot(file):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'ionbot' ... ", end = '', flush = True)
    df = pd.read_csv(file)

    df = df[["uniprot_id", "modification", "position"]]

    # get rid of prefixes in protein ids
    df["uniprot_id"] = df["uniprot_id"].str.replace("sp|", "", regex=False).str.replace("tr|", "", regex=False)

    # fix ionbot bug where C-term modifications are not correctly placed in the protein
    df['position'] = df.apply(lambda x: x['position'] - 1 if 'C-term' in x['modification'] else x['position'], axis=1)

    df["modification_unimod_name"] = df["modification"].apply(
        lambda x: x.split("]")[1].split("[")[0].lower()
    )

    df["modification_unimod_id"] = df["modification"].apply(
        lambda x: int(x.split("]")[0][1:])
    )

    df["modified_residue"] = df["modification"].apply(
        lambda x: x.split("[")[2].split("]")[0]
    )

    #get all unique unimod id - amino acid combinations
    unique_mods = df[["modification_unimod_id", "modified_residue"]].drop_duplicates()

    #look unimod classification up in the db, store them in column
    unique_mods["classification"] = unique_mods.apply(
        lambda x: classification_from_id(
            x["modification_unimod_id"], x["modified_residue"], UNIMOD_MAPPER
        ),
        axis=1,
    )

    #merge
    df = df.merge(unique_mods, on=["modification_unimod_id", "modified_residue"], how="left")

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    #get all unique unimod ids
    unique_mass_shifts = df[["modification_unimod_id"]].drop_duplicates()

    #look mass shifts up in the db, store them in column
    unique_mass_shifts["mass_shift"] = unique_mass_shifts["modification_unimod_id"].apply(
        lambda x: mass_shift_from_id(x, UNIMOD_MAPPER)
    )

    #look unimod names up in the db, store them in column (ionbot renames some modifications)
    unique_mass_shifts["modification_unimod_name_mapped"] = unique_mass_shifts["modification_unimod_id"].apply(
        lambda x: map_id_to_unimod_name(x)
    )

    #merge
    df = df.merge(unique_mass_shifts, on="modification_unimod_id", how="left")

    # use mapped name if available, else use original name
    df["modification_unimod_name"] = df.apply(
        lambda row: row["modification_unimod_name"] if pd.isna(row["modification_unimod_name_mapped"]) else row["modification_unimod_name_mapped"],
        axis=1)

    # get mismatch candidates
    df["candidates"] = df.apply(
        lambda row: get_candidates_from_mass_shift(row["mass_shift"]),
        axis=1
    )

    df["display_name"] = df["modification_unimod_name"]

    df.drop(columns=["modification", "modified_residue"], inplace=True)

    return df


def read_mod_csv(file):
    """
    read plain csv.
    minimum columns: uniprot_id, position, modification_unimod_id
    optional: classification, mass_shift, annotation
    if classification and mass_shift are not provided, they will be queried from unimod
    """
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'csv' ... ", end = '', flush = True)
    df = pd.read_csv(file)
    if "classification" not in [x.lower() for x in df.columns]:
        if "modified_residue" not in [x.lower() for x in df.columns]:
            fasta = query_uniprot_for_fasta(df["uniprot_id"].tolist())
            df["modified_residue"] = df.apply(
                lambda x: get_amino_acid_from_protein(
                    x["position"], x["uniprot_id"], fasta
                ),
                axis=1,
            )

        df["classification"] = df.apply(
            lambda x: classification_from_id(
                x["modification_unimod_id"], x["modified_residue"], UNIMOD_MAPPER
            ),
            axis=1,
        )

        df.drop(columns=["modified_residue"], inplace=True)

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    if "mass_shift" not in [x.lower() for x in df.columns]:
        df["mass_shift"] = df["modification_unimod_id"].apply(
            lambda x: mass_shift_from_id(x, UNIMOD_MAPPER)
        )

    if "display_name" not in [x.lower() for x in df.columns]:
        df["display_name"] = df["modification_unimod_name"]

    return df


def read_sage(file):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'sage' ... ", end = '', flush = True)
    # write sage to temporary file and then give filepath to psm utils
    with NamedTemporaryFile(mode="wt", delete=False) as f:
        f.write(file.getvalue())
        f.close()

    psm_list = read_file(f.name, filetype="sage")  # file is a StringIO object!
    Path(f.name).unlink()

    psm_list = psm_list[psm_list["qvalue"] <= 0.01]  # filter at 1% FDR
    psm_list = [x for x in psm_list if not x.is_decoy]  # remove decoys
    psm_list = [
        x for x in psm_list if x.peptidoform.is_modified
    ]  # filter for modified peptidoforms
    psm_list = [
        x for x in psm_list if len(x.protein_list) == 1
    ]  # filter nonunique peptidoforms

    df = PSMList_to_mod_df(
        psm_list
    )  # [protein, mass_shift, modification_unimod_id, modification_unimod_name, display_name, position]

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    return df


def read_maxquant(file):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'maxquant' ... ", end = '', flush = True)
    # write msms to temporary file and then give filepath to psm utils
    with NamedTemporaryFile(mode="wt", delete=False) as f:
        f.write(file.getvalue())
        f.close()

    psm_list = read_file(f.name, filetype="msms")  # file is a StringIO object!
    Path(f.name).unlink()

    maxquant_mapper = pd.read_csv(BASEPATH + "/static/resources/map_mq_file.csv", index_col=0)[
        "value"
    ].to_dict()
    psm_list.rename_modifications(maxquant_mapper)

    psm_list = [
        x for x in psm_list if x.peptidoform.is_modified
    ]  # filter for modified peptidoforms
    psm_list = [
        x for x in psm_list if len(x.protein_list) == 1
    ]  # filter nonunique peptidoforms

    df = PSMList_to_mod_df(psm_list)

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    return df


def read_spectronaut(file):
    # reads spectronaut PTM Site report file
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'spectronaut' ... ", end = '', flush = True)

    df = pd.read_csv(file, usecols=["PTM.ProteinId", "PTM.SiteLocation", "PTM.ModificationTitle", "PTM.SiteAA"], sep="\t")
    df = df.dropna(subset=["PTM.SiteLocation"])
    df["PTM.SiteLocation"] = df["PTM.SiteLocation"].astype(int)
    df.rename(columns={"PTM.ProteinId": "uniprot_id", "PTM.SiteLocation": "position", "PTM.ModificationTitle": "modification_unimod_name", "PTM.SiteAA": "modified_residue"}, inplace=True)
    df["modification_unimod_name"] = df["modification_unimod_name"].apply(lambda x: x.split(" (")[0])

    unique_mods = df[["modification_unimod_name", "modified_residue"]].drop_duplicates()

    unique_mods["modification_unimod_id"] = unique_mods["modification_unimod_name"].apply(
        lambda x: map_modification_string_to_unimod_id(x)
    )
    unique_mods["classification"] = unique_mods.apply(
        lambda x: classification_from_id(x["modification_unimod_id"], x["modified_residue"], UNIMOD_MAPPER),
        axis=1
    )

    unique_mods["mass_shift"] = unique_mods["modification_unimod_id"].apply(
        lambda x: mass_shift_from_id(int(x), UNIMOD_MAPPER)
    )

    df = df.merge(unique_mods, on=["modification_unimod_name", "modified_residue"], how="left")

    df["display_name"] = df["modification_unimod_name"]

    df = df[
        [
            "uniprot_id",
            "position",
            "modification_unimod_name",
            "modification_unimod_id",
            "classification",
            "mass_shift",
            "display_name",
        ]
    ]

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    return df


def read_any(file):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'any' ... ", end = '', flush = True)
    # write to temporary file and then give filepath to psm utils

    with NamedTemporaryFile(mode="wt", delete=False) as f:
        f.write(file.getvalue())
        f.close()

    psm_list = read_file(f.name, filetype="infer")  # file is a StringIO object!
    Path(f.name).unlink()

    psm_list = psm_list[psm_list["qvalue"] <= 0.01]  # filter at 1% FDR
    psm_list = [x for x in psm_list if not x.is_decoy]  # remove decoys
    psm_list = [
        x for x in psm_list if x.peptidoform.is_modified
    ]  # filter for modified peptidoforms
    psm_list = [
        x for x in psm_list if len(x.protein_list) == 1
    ]  # filter nonunique peptidoforms

    df = PSMList_to_mod_df(
        psm_list
    )  # [protein, mass_shift, modification_unimod_id, modification_unimod_name, display_name, position]

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    return df


def read_mzid(file):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tProcess data of type 'mzid' ... ", end = '', flush = True)
    # write mzid to temporary file and then give filepath to psm utils
    with NamedTemporaryFile(mode="wt", delete=False) as f:
        f.write(file.getvalue())
        f.close()

    psm_list = read_file(f.name, filetype="mzid")  # file is a StringIO object!
    Path(f.name).unlink()

    psm_list = [x for x in psm_list if not x.is_decoy]  # remove decoys

    psm_list = [
        x for x in psm_list if x.peptidoform.is_modified
    ]  # filter for modified peptidoforms

    psm_list = [
        x for x in psm_list if len(x.protein_list) == 1
    ]  # filter nonunique peptidoforms

    df = PSMList_to_mod_df(
        psm_list
    )  # [protein, mass_shift, modification_unimod_id, modification_unimod_name, display_name, position]

    #remove unwanted modifications based on unimod class assignment
    global EXCLUDE_CLASSES
    df = df[~df["classification"].isin(EXCLUDE_CLASSES)]

    return df


def construct_modifications_entry(row):
    return {
        "modification_unimod_name": row["modification_unimod_name"],
        "modification_classification": row["classification"]
        if "classification" in row
        else "null",
        "modification_unimod_id": row["modification_unimod_id"],
        "display_name": row["display_name"] if "display_name" in row else "null",
        "mass_shift": row["mass_shift"] if "mass_shift" in row else "null",
        "candidates": row["candidates"] if "candidates" in row else "null",
    }


def parse_df_to_json_schema(dataframe):
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tConvert data into schema ... ", end = '', flush = True)
    # Init. dictionary to store results in JSON schema.
    protein_dict = {"proteins": {}, "meta_data": {}}
    # For each row
    for _, row in dataframe.iterrows():
        # Add entire entry, if protein does not exist yet.
        if row["uniprot_id"] not in protein_dict["proteins"].keys():
            entry = {
                row["uniprot_id"]: {
                    "positions": {
                        row["position"]: {
                            "modifications": [construct_modifications_entry(row)]
                        }
                    },
                }
            }
            protein_dict["proteins"].update(entry)
        # Add position and modification identifier, name, class, display name and annotations, if protein already exists but not the position.
        elif (
            row["position"]
            not in protein_dict["proteins"][row["uniprot_id"]]["positions"].keys()
        ):
            protein_dict["proteins"][row["uniprot_id"]]["positions"][
                row["position"]
            ] = {"modifications": [construct_modifications_entry(row)]}

        # Add modification identifier, name, class and annotations, if the position in the protein already exists.
        elif (
            row["position"]
            in protein_dict["proteins"][row["uniprot_id"]]["positions"].keys()
        ):
            protein_dict["proteins"][row["uniprot_id"]]["positions"][row["position"]][
                "modifications"
            ].append(construct_modifications_entry(row))
    if DEBUG :
        print( "DONE" )
    return protein_dict


def read_userinput(file, flag):
    # psm-utils based parsing (sage, msms) might not work for any other encoding than UTF-8, see psm-utils issue #68
    if flag == "ionbot":
        df = read_ionbot(file)
    elif flag == "msfragger":
        df = read_msfragger(file)
    elif flag == "csv":
        df = read_mod_csv(file)
    elif flag == "sage":
        df = read_sage(file)
    elif flag == "msms":
        df = read_maxquant(file)
    elif flag == "spectronaut":
        df = read_spectronaut(file)
    elif flag == "mzid":
        df = read_mzid(file)
    elif flag == "infer":
        df = read_any(file)
    df = df.fillna("null")
    if DEBUG :
        print( "DONE" )
    return parse_df_to_json_schema(df.drop_duplicates())


def parse_user_input(user_file, user_flag, mass_shift_tolerance = 0.001, exclude_classes = []):
    global TOLERANCE
    TOLERANCE = mass_shift_tolerance
    global EXCLUDE_CLASSES
    EXCLUDE_CLASSES = exclude_classes
    try:
        json = read_userinput(user_file, user_flag)
        json["meta_data"]["mass_shift_tolerance"] = mass_shift_tolerance
        json["meta_data"]["exclude_classes"] = exclude_classes
    except TypeError:
        raise TypeError(
            "Your file could not be parsed. Have you selected the right format?"
        )
    return json


def _brotli_decompress(content: str) -> str:
    return brotli.decompress(base64.urlsafe_b64decode(content)).decode()


def _brotly_compress(content: str) -> str:
    return base64.urlsafe_b64encode(brotli.compress(content.encode())).decode()