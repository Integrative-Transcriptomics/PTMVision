from Bio.PDB import PDBParser
from Bio import SeqIO
from io import StringIO
from unimod_mapper import UnimodMapper
import pandas as pd
import numpy as np
import scipy as sc
from pathlib import Path
import re
import urllib.request
import requests
import brotli
import base64
from psm_utils.io import read_file
from tempfile import NamedTemporaryFile
from pyteomics.mass.unimod import Unimod
import json

MAPPER = UnimodMapper()
PDBPARSER = PDBParser(PERMISSIVE=False)
BASEPATH = "/app/ptmvision"

"""
TODO: Refactor into reader classes, one for each file format.
"""

def check_fasta_coverage(uniprot_ids, fasta_dict):
    # check how many proteins were found in uniprot, let the user decide if they want to provide their own fasta
    found = len(fasta_dict)
    total = len(uniprot_ids)
    not_found = total - found
    if not_found > 50:
        raise Exception(
            "Not able to retrieve protein sequences for {} protein IDs. Please provide the fasta used for database search.".format(
                not_found
            )
        )


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
            contacts.setdefault(int(pairwise_contact[0]), []).append(
                int(pairwise_contact[1])
            )
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


def get_ptm_position_in_protein(mod, protein_start):
    if protein_start == -1:
        return -1
    position_ptm_in_peptide = int(re.search(".+?(?=\D)", mod)[0])
    position_ptm_in_protein = protein_start + position_ptm_in_peptide - 1
    return position_ptm_in_protein


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


def map_mass_to_unimod_msfragger(mods):
    modstr = ""
    for mod in mods.split(","):
        pos_aa = mod.split("(")[0]

        mass = re.search(r"\(.*?\)", mod)[0]
        mass = float(mass.replace("(", "").replace(")", ""))
        mapped = MAPPER.mass_to_names(
            mass, decimals=4
        )  # for unimod ID: mapper.mass_to_ids(mass, decimals = 4)

        mapped_names = " or ".join(list(mapped))
        if len(list(mapped)) > 0:
            modstr += pos_aa + "[" + mapped_names + "]"

    return modstr


def map_mass_to_unimod_name(mass_shift):
    mapped = MAPPER.mass_to_names(mass_shift, decimals=4)
    if len(mapped) > 1:  # mass shift could not be uniquely mapped to unimod
        return " or ".join(list(mapped))
    elif len(mapped) == 0:  # mass shift couldnt bet mapped at all
        return None
    else:
        return mapped[0]


def map_unimod_name_to_mass(unimod_name):
    mapped = MAPPER.name_to_mass(unimod_name)
    if len(mapped) == 0:
        return None
    else:
        return mapped[0]


def map_id_to_unimod_name(id):
    mapped = list(set(MAPPER.id_to_name(id)))
    if len(mapped) > 1:  # name could not be uniquely mapped to unimod
        return " or ".join(list(mapped))
    elif len(mapped) == 0:  # mass shift couldnt bet mapped at all
        return None
    else:
        return mapped[0]


def map_id_to_mass(id):
    mapped = MAPPER.id_to_mass(id)
    if len(mapped) > 1:
        mass_shifts = [str(x) for x in list(set(mapped))]
        return " or ".join(mass_shifts)
    elif len(mapped) == 0:
        return None
    else:
        return mapped[0]


def map_mass_to_unimod_id(mass_shift):
    mapped = MAPPER.mass_to_ids(mass_shift, decimals=4)
    if len(mapped) > 1:  # mass shift could not be uniquely mapped to unimod
        return " or ".join(list(mapped))
    elif len(mapped) == 0:  # mass shift couldnt bet mapped at all
        return None
    else:
        return mapped[0]


def map_unimod_description_to_unimod_id(mod):
    query_string = "`Description` == '{}'".format(mod)
    matches = list(set(MAPPER.query(query_string)["Accession"].tolist()))
    if len(matches) > 1:  # mass shift could not be uniquely mapped to unimod
        return " or ".join(list(matches))
    elif len(matches) == 0:  # mass shift couldn't bet mapped at all
        return None
    else:
        return matches[0]


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
    try:
        mod = unimod_db.get(int(unimod_id))
    except KeyError:  # deprecated unimod ID :(
        return 
    for specification in mod.specificities:
        if specification.amino_acid == amino_acid:
            classification = specification.classification.classification
            return classification

    return 


def get_classification(unimod_id, residue, unimod_db):
    if unimod_id is None:
        return
    elif " or " in unimod_id:
        return
    else:
        classification = classification_from_id(unimod_id, residue, unimod_db)
        return classification


def mass_shift_from_id(unimod_id, unimod_db):
    """
    Get mass shift from unimod ID
    """
    try:
        mod = unimod_db.get(unimod_id)
    except KeyError:  # deprecated unimod ID :(
        return None
    return mod.monoisotopic_mass


def id_from_name(name):
    ids = []
    for mod_name in name.split(" or "):
        id = MAPPER.name_to_id(mod_name)[0]
        ids.append(id)
    return ",".join(ids)


def from_psm_to_protein(df):
    # for each psm in msfragger, map PTM position onto protein
    uniprot_ids = []
    modifications = []
    positions = []
    for index, psm in df.iterrows():
        protein_start = psm["Protein Start"]
        uniprot_id = psm["Protein ID"]

        for mod in psm["Assigned Modifications parsed"].split("]"):
            if mod == "":
                continue
            modification = mod.split("[")[1]
            position = get_ptm_position_in_protein(mod, protein_start)

            modifications.append(modification)
            positions.append(position)
            uniprot_ids.append(uniprot_id)

    zipped = list(zip(uniprot_ids, positions, modifications))
    protein_df = pd.DataFrame(
        zipped, columns=["uniprot_id", "position", "modification_unimod_name"]
    )
    protein_df["modification_unimod_id"] = protein_df["modification_unimod_name"].apply(
        lambda x: id_from_name(x)
    )

    return protein_df


def extract_mods_from_proforma(peptidoform):
    """
    ENYC[+57.0215]NNVMM[+15.994915]K/2	-> [(3, 57.0215), (8, 15.993915)] (zero indexed!), try to map to unimod names
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


def parse_protein_string_old(x):
    """
    Try to parse uniprot ID from protein_list entry in PSMList
    """
    uniprot_accession_pattern = re.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )
    accession_candidates = [x]
    if "|" in x:
        if x.count("|") == 2:
            accession_candidates.append(x.split("|")[1])
        else:
            return accession_candidates.append(x.split("|")[-1])
    for x in accession_candidates:
        if uniprot_accession_pattern.match(x):
            return x  # probably accession number
    return None


def parse_protein_string(x):
    """
    Try to parse uniprot accession ID from protein string
    """
    uniprot_accession_pattern = re.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    )
    accession_candidates = [x] + x.split("|")
    for x in accession_candidates:
        if uniprot_accession_pattern.match(x):
            return x  # probably accession number
    return None


def msfragger_to_unimod_id(x):
    """
    map mass shift or description to unimod name
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
    if row["modification_unimod_name"]:
        return row["modification_unimod_name"]
    else:
        return row["modification"][1]


def get_mass_shift(mod):
    # gets with either name ("Oxidation") or string ("+57.") with mass shift
    if mod[0] in ["+", "-"]:
        mass_shift = float(mod)
    else:
        # try to get mass over unimod mapper
        mass_shift = map_unimod_name_to_mass(mod)
    return mass_shift


def PSMList_to_mod_df(psm_list):
    """
    parse list of PSMs and return pandas dataframe with [protein, position, modification_unimod_id, modification_unimod_name]
    """
    peptidoforms = []
    proteins = []
    peptide_sequences = []
    for psm in psm_list:
        peptidoforms.append(psm["peptidoform"].proforma)
        proteins.append(psm.protein_list[0])
        peptide_sequences.append(psm["peptidoform"].sequence)

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

    # query uniprot for protein sequence
    fasta = query_uniprot_for_fasta(df["uniprot_id"].tolist())

    df["protein_sequence"] = df["uniprot_id"].apply(
        lambda x: get_protein_sequence(fasta, x)
    )

    df = df[df["protein_sequence"] != ""]

    df["mass_shift"] = df["modification"].apply(lambda x: get_mass_shift(x[1]))
    df["mod_position_in_peptide"] = df["modification"].apply(lambda x: x[0])

    # TODO: in case of already annotated mass shift, parse unimod name and id from modification col directly
    df["modification_unimod_name"] = df["mass_shift"].apply(
        lambda x: map_mass_to_unimod_name(x)
    )
    df["modification_unimod_id"] = df["mass_shift"].apply(
        lambda x: map_mass_to_unimod_id(x)
    )
    df["display_name"] = df.apply(
        lambda x: get_display_name(x), axis=1
    )  # Unimod Name if available, otherwise "unannotated mass shift: (mass shift)"
    df["peptide_position"] = df.apply(
        lambda x: get_peptide_position_in_protein(x["peptide"], x["protein_sequence"]),
        axis=1,
    )
    df["position"] = df["mod_position_in_peptide"] + df["peptide_position"]

    unimod_db = Unimod()
    df["classification"] = df.apply(
        lambda x: get_classification(
            x["modification_unimod_id"],
            x["peptide"][x["mod_position_in_peptide"]],
            unimod_db,
        ),
        axis=1,
    )

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
        mod_list.append((pos, "Unannotated mass-shift " +str(mass)))

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
    parse modifications from MSFragger PSM
    variable mods: Assigned Modification ("7M(15.9949)")  N-term(42.0106)
    new mods: Observed Modiciation (Mod1: Methylation (PeakApex: 14.0144, Theoretical: 14.0157), localization in MSFragger Localization (TVkEAEEAAK)
    combine into one list of tuples => ((7,15.9949), ())
    """
    mods = []
    if not pd.isna(msfragger_psm_row["Assigned Modifications"]):
        var_mods = parse_assigned_mods_msfragger(
            msfragger_psm_row["Assigned Modifications"]
        )
        if var_mods:
            mods = mods + var_mods
    if not pd.isna(msfragger_psm_row["MSFragger Localization"]):
        new_mods = parse_observed_mods_msfragger(
            msfragger_psm_row["MSFragger Localization"],
            msfragger_psm_row["Observed Modifications"],
        )
        if new_mods:
            mods = mods + new_mods

    return mods


def read_msfragger(file):
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
    pept_mods = pept_mods[pept_mods["Observed Modifications"].str.count(";") == 0]

    # filter out combinations of modifications
    pept_mods = pept_mods[pept_mods["Observed Modifications"].str.count("Mod2") == 0]

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

    # Protein Start is 1 indexed, modification position is 0 indexed, so cancels out the -1
    pept_mods = pept_mods.dropna(subset=["modification"])
    pept_mods["modification_unimod_id"] = pept_mods["modification"].apply(
        lambda x: msfragger_to_unimod_id(x)
    )
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
    unimod_db = Unimod()
    pept_mods["residue"] = pept_mods.apply(
        lambda x: x["Peptide"][int(x["modification"][0])], axis=1
    )
    pept_mods["classification"] = pept_mods.apply(
        lambda x: get_classification(
            x["modification_unimod_id"], x["residue"], unimod_db
        ),
        axis=1,
    )

    # get display name (unimod name if available, otherwise "unannotated mass shift: (mass shift)"
    pept_mods["display_name"] = pept_mods.apply(lambda x: get_display_name(x), axis=1)

    # pept_mods = pept_mods[
    #    [
    #        "uniprot_id",
    #        "position",
    #        "modification_unimod_name",
    #        "modification_unimod_id",
    #        "classification",
    #        "mass_shift",
    #        "display_name",
    #    ]
    # ]

    return pept_mods


def read_ionbot(file):
    df = pd.read_csv(file)
    unimod_db = Unimod()

    df = df[["uniprot_id", "modification", "position"]]
    df["modification_unimod_name"] = df["modification"].apply(
        lambda x: x.split("]")[1].split("[")[0].lower()
    )
    df["modification_unimod_id"] = df["modification"].apply(
        lambda x: int(x.split("]")[0][1:])
    )
    df["modified_residue"] = df["modification"].apply(
        lambda x: x.split("[")[2].split("]")[0]
    )
    df["classification"] = df.apply(
        lambda x: classification_from_id(
            x["modification_unimod_id"], x["modified_residue"], unimod_db
        ),
        axis=1,
    )
    df["mass_shift"] = df["modification_unimod_id"].apply(
        lambda x: mass_shift_from_id(x, unimod_db)
    )
    df["display_name"] = df["modification_unimod_name"]

    df.drop(columns=["modification", "modified_residue"], inplace=True)
    
    df["uniprot_id"] = df["uniprot_id"].str.replace("sp|", "", regex=False).str.replace("tr|", "", regex=False)

    return df


def read_mod_csv(file):
    """
    read plain csv.
    minimum columns: uniprot_id, position, modification_unimod_id
    optional: classification, mass_shift, annotation
    if classification and mass_shift are not provided, they will be queried from unimod
    """
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

        unimod_db = Unimod()
        df["classification"] = df.apply(
            lambda x: classification_from_id(
                x["modification_unimod_id"], x["modified_residue"], unimod_db
            ),
            axis=1,
        )
        df.drop(columns=["modified_residue"], inplace=True)
    if "mass_shift" not in [x.lower() for x in df.columns]:
        unimod_db = Unimod()
        df["mass_shift"] = df["modification_unimod_id"].apply(
            lambda x: mass_shift_from_id(x, unimod_db)
        )
    if "display_name" not in [x.lower() for x in df.columns]:
        df["display_name"] = df["modification_unimod_name"]

    return df


def read_sage(file):
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

    return df


def read_maxquant(file):
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

    return df


def read_any(file):
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
    }


def parse_df_to_json_schema(dataframe):
    # Init. dictionary to store results in JSON schema.
    protein_dict = {"proteins": {}, "meta_data": {}}
    # Internal function to add entries.
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
    elif flag == "infer":
        df = read_any(file)

    df = df.fillna("null")
    
    return parse_df_to_json_schema(df.drop_duplicates())


def parse_user_input(user_file, user_flag):
    try:
        json = read_userinput(user_file, user_flag)
    except TypeError:
        raise TypeError(
            "Your file could not be parsed. Have you selected the right format?"
        )
    return json


def _brotli_decompress(content: str) -> str:
    return brotli.decompress(base64.urlsafe_b64decode(content)).decode()


def _brotly_compress(content: str) -> str:
    return base64.urlsafe_b64encode(brotli.compress(content.encode())).decode()


if __name__ == "__main__":
    # json = parse_user_input("example_data/sage_results_v0131.tsv", "sage")
    json = parse_user_input(
        "example_data/msfragger_opensearch_ptmshepherdv206_msfraggerv40_psm.tsv",
        "msfragger",
    )
    # json = parse_user_input("example_data/msms_maxquant_v24130.txt", "msms")
    # json = parse_user_input("example_data/sulfolobus_acidocaldarius_UP000060043_with_crap.sage_default_config.results.sage.tsv", "sage")
    # print(json)
