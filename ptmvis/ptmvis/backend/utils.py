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

mapper = UnimodMapper()
pdbparser = PDBParser(PERMISSIVE=False)


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
        structure = pdbparser.get_structure(uniprot_id, StringIO(pdb_text))
    except:
        raise Exception("Alphafold structure not available.")
    return structure, pdb_text


def parse_structure(structure_string):
    structure = None
    # with warnings.filterwarnings("error"):
    try:
        structure = pdbparser.get_structure("custom_pdb", StringIO(structure_string))
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


def sage_peptide_to_msfragger_mod(sage_peptide):
    modified_sites = re.finditer("\[", sage_peptide)
    indices = [index.start() for index in modified_sites]
    assigned_modifications = []

    for index in indices:
        modified_aa = sage_peptide[index - 1]
        weight = sage_peptide[index + 1 :].split("]")[0].replace("+", "")

        # get index of modified peptide (not same as index bc modifications are in the string as well)
        prefix = sage_peptide[:index]
        prefix = (
            re.sub(r"\[(.*?)\]", "", prefix)
            .replace("[", "")
            .replace("]", "")
            .replace("-", "")
        )
        mod_index = len(prefix)
        assigned_modifications.append(str(mod_index) + modified_aa + "(" + weight + ")")

    return ",".join(assigned_modifications)


def map_mass_to_unimod_msfragger(mods):
    modstr = ""
    for mod in mods.split(","):
        pos_aa = mod.split("(")[0]

        mass = re.search(r"\(.*?\)", mod)[0]
        mass = float(mass.replace("(", "").replace(")", ""))
        mapped = mapper.mass_to_names(
            mass, decimals=4
        )  # for unimod ID: mapper.mass_to_ids(mass, decimals = 4)

        mapped_names = " or ".join(list(mapped))
        if len(list(mapped)) > 0:
            modstr += pos_aa + "[" + mapped_names + "]"

    return modstr


def map_mass_to_unimod_name(mass_shift):
    mapped = mapper.mass_to_names(mass_shift, decimals=4)
    if len(mapped) > 1:  # mass shift could not be uniquely mapped to unimod
        return " or ".join(list(mapped))
    elif len(mapped) == 0:  # mass shift couldnt bet mapped at all
        return None
    else:
        return mapped[0]


def map_unimod_name_to_mass(unimod_name):
    mapped = mapper.name_to_mass(unimod_name)
    if len(mapped) == 0:
        return None
    else:
        return mapped[0]


def map_mass_to_unimod_id(mass_shift):
    mapped = mapper.mass_to_ids(mass_shift, decimals=4)
    if len(mapped) > 1:  # mass shift could not be uniquely mapped to unimod
        return " or ".join(list(mapped))
    elif len(mapped) == 0:  # mass shift couldnt bet mapped at all
        return None
    else:
        return mapped[0]


def id_from_name(name):
    ids = []
    for mod_name in name.split(" or "):
        id = mapper.name_to_id(mod_name)[0]
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
    ENYC[+57.0215]NNVMM[+15.994915]K/2	-> ((3, 57.0215), (8, 15.993915)) (zero indexed!), try to map to unimod names
    AADM[Oxidation]TGADIEAMTR/2 -> ((3, Oxidation))
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


def get_display_name(row):
    if row["modification_unimod_name"]:
        return row["modification_unimod_name"]
    else:
        return "unannotated mass shift {}".format(row["mass_shift"])


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

    df = df[
        [
            "uniprot_id",
            "mass_shift",
            "modification_unimod_name",
            "modification_unimod_id",
            "display_name",
            "position",
        ]
    ]
    df = df.drop_duplicates()

    return df


def read_msfragger(file):
    # read file and pick relevant rows and columns
    print("Reading MSFragger file...")
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

    pept_mods = pept_mods[
        ~(pept_mods["Assigned Modifications"].isna())
        | (pept_mods["MSFragger Localization"].isna())
    ]

    # try to assign modifications to mass shifts
    pept_mods["Assigned Modifications parsed"] = pept_mods[
        "Assigned Modifications"
    ].apply(lambda x: map_mass_to_unimod_msfragger(x))

    # rewrite to protein level modifications
    print("Mapping modification sites onto proteins...")
    prot_mods = from_psm_to_protein(pept_mods)

    return prot_mods


def read_ionbot(file):
    df = pd.read_csv(file)
    df = df[["uniprot_id", "unexpected_modification", "position"]]
    df["modification_unimod_name"] = df["unexpected_modification"].apply(
        lambda x: x.split("]")[1].split("[")[0].lower()
    )
    df["modification_unimod_id"] = df["unexpected_modification"].apply(
        lambda x: x.split("]")[0][1:]
    )
    df.drop(columns=["unexpected_modification"], inplace=True)
    return df


def read_mod_csv(file):
    df = pd.read_csv(file)
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

    maxquant_mapper = pd.read_csv("static/resources/map_mq_file.csv", index_col=0)[
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


def parse_df_to_json_schema(dataframe):
    # Init. dictionary to store results in JSON schema.
    protein_dict = {"proteins": {}, "meta_data": {}}
    # Internal function to add entries.
    construct_modifications_entry = lambda row: {
        "modification_unimod_name": row["modification_unimod_name"],
        "modification_classification": row["classification"]
        if "classification" in row
        else "null",
        "modification_unimod_id": row["modification_unimod_id"],
    }
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
                    "pdb_structure": "null",
                }
            }
            protein_dict["proteins"].update(entry)
        # Add position and modification identifier, name, class and annotations, if protein already exists but not the position.
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
    return parse_df_to_json_schema(df)


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
    json = parse_user_input("example_data/msms_maxquant_v24130.txt", "msms")
    # json = parse_user_input("example_data/sulfolobus_acidocaldarius_UP000060043_with_crap.sage_default_config.results.sage.tsv", "sage")
    print(json)
