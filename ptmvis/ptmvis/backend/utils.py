import re, warnings, urllib.request
import pandas as pd
import numpy as np
import scipy as sc
from Bio.PDB import PDBParser
from Bio import SeqIO
from io import StringIO
from unimod_mapper import UnimodMapper
import requests
import json

mapper = UnimodMapper()
pdbparser = PDBParser(PERMISSIVE=False)


def get_distance_matrix(structure):
    residue_count = 0
    for residue in enumerate(structure.get_residues()):
        residue_count += 1
    calpha_coords = np.zeros(shape=(residue_count, 3))

    for i, residue in enumerate(structure.get_residues()):
        for atom in residue:
            if atom.get_name() == "CA":
                calpha_coords[i] = atom.get_coord()
                break
    residue_distances = sc.spatial.distance_matrix(calpha_coords, calpha_coords)
    return residue_distances


def get_contacts(dist_matrix, threshold):
    contacts = {}
    for pairwise_contact in np.argwhere(dist_matrix < threshold):
        if pairwise_contact[0] != pairwise_contact[1]:
            contacts.setdefault(int(pairwise_contact[0]), []).append(
                int(pairwise_contact[1])
            )
    return contacts


def get_structure(unimod_id):
    url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb".format(unimod_id)
    structure = None
    try:
        response = urllib.request.urlopen(url)
        data = response.read()
        text = data.decode("utf-8")
        structure = pdbparser.get_structure(unimod_id, StringIO(text))
    except:
        raise Exception("Alphafold structure not available.")
    return structure


def parse_structure(structure_string):
    structure = None
    with warnings.filterwarnings("error"):
        try:
            structure = pdbparser.get_structure(
                "custom_pdb", StringIO(structure_string)
            )
        except Exception as e:
            raise Exception("Failed to parse provided .pdb structure: " + str(e))
        else:
            return structure


def get_fasta(uniprot_ids):

    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28"
    url_separator = "%20OR%20"
    id_list = ["%28id%3A{}".format(uniprot_id + "%29") for uniprot_id in uniprot_ids]

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
        fasta_dict[rec.id.split("|")[-1]] = str(rec.seq)
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


def map_mass_to_unimod(mods):
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


def id_from_name(name):
    ids = []
    for mod_name in name.split(" or "):
        id = mapper.name_to_id(mod_name)[0]
        ids.append(id)
    return ",".join(ids)


def from_psm_to_protein(df):
    # for each psm in msfragger, get protein sequence via uniprot ID and map PTM position onto it
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
            "Protein Start",
        ]
    ]
    pept_mods = pept_mods.dropna(subset=["Modified Peptide"])
    pept_mods = pept_mods.drop_duplicates()

    # try to assign modifications to mass shifts
    print("Assigning modifications to mass shifts...")
    pept_mods["Assigned Modifications parsed"] = pept_mods[
        "Assigned Modifications"
    ].apply(lambda x: map_mass_to_unimod(x))

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
    json = parse_df_to_json_schema(df)
    return df


def read_sage_csv(file):
    df = pd.read_table(file)
    df.columns = [
        x.lower() for x in df.columns
    ]  # sometimes its upper, sometimes its lower case (version specific?)

    df = df[["peptide", "proteins", "label"]]

    # drop decoy matches
    df = df.loc[df["label"] == 1]
    df = df.drop(columns=["label"])

    # drop PSMs without modification
    df = df.loc[df["peptide"].str.contains("\[")]

    # one protein per line
    df = (
        df.set_index(["peptide"])
        .apply(lambda x: x.str.split(";").explode())
        .reset_index()
    )

    # parse uniprot Accession IDs (but use IDs for uniprot mapping)
    df["Protein ID"] = df["proteins"].apply(lambda x: x.split("|")[-2])

    # get protein sequences (first: try uniprot, if it fails, let user upload fasta)
    uniprot_ids = df["proteins"].unique().tolist()
    uniprot_ids = [x.split("|")[-1] for x in uniprot_ids]
    fasta_dict = get_fasta(uniprot_ids)

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

    # add protein sequences to df
    df["protein_sequence"] = df["proteins"].apply(
        lambda x: get_protein_sequence(fasta_dict, x.split("|")[-1])
    )
    df = df.loc[df["protein_sequence"] != ""]

    # find Protein Start index of peptide
    df["Protein Start"] = df.apply(
        lambda x: get_peptide_position_in_protein(x["peptide"], x["protein_sequence"]),
        axis=1,
    )

    # convert to MSFragger Assigned Modifications Format (10S(79.9663), 8M(15.9949)) so we can reuse functions
    df["Assigned Modifications"] = df["peptide"].apply(
        lambda x: sage_peptide_to_msfragger_mod(x)
    )

    # try to assign modifications to mass shifts
    print("Assigning modifications to mass shifts...")
    df["Assigned Modifications parsed"] = df["Assigned Modifications"].apply(
        lambda x: map_mass_to_unimod(x)
    )

    # rewrite to protein level modifications
    print("Mapping modification sites onto proteins...")
    df = from_psm_to_protein(df)

    return df


def read_file(file, flag):
    if flag == "ionbot":
        df = read_ionbot(file)
    elif flag == "msfragger":
        df = read_msfragger(file)
    elif flag == "csv":
        df = read_mod_csv(file)
    elif flag == "sage":
        df = read_sage_csv(file)
    return parse_df_to_json_schema(df)


def parse_df_to_json_schema(dataframe):
    protein_dict = {"proteins": {}, "meta_data": ()}

    for i, row in dataframe.iterrows():
        # if protein is not there, add entire entry
        if row["uniprot_id"] not in protein_dict["proteins"].keys():
            entry = {
                row["uniprot_id"]: {
                    "positions": {
                        row["position"]: {
                            "modifications": [
                                {
                                    "modification_unimod_name": row[
                                        "modification_unimod_name"
                                    ],
                                    "modification_classification": row[
                                        "classification"
                                    ],
                                    "modification_unimod_id": row[
                                        "modification_unimod_id"
                                    ],
                                }
                            ]
                        }
                    },
                    "pdb_structure": "NaN",
                }
            }
            protein_dict["proteins"].update(entry)

        # if position is not there, but protein is, add position and mod id/name/class/annotations
        elif (
            row["position"]
            not in protein_dict["proteins"][row["uniprot_id"]]["positions"].keys()
        ):
            protein_dict["proteins"][row["uniprot_id"]]["positions"][
                row["position"]
            ] = {
                "modifications": [
                    {
                        "modification_unimod_name": row["modification_unimod_name"],
                        "modification_classification": row["classification"],
                        "modification_unimod_id": row["modification_unimod_id"],
                    }
                ]
            }

        # if position in this protein is already in dict, add mod id/name/class/annotations
        elif (
            row["position"]
            in protein_dict["proteins"][row["uniprot_id"]]["positions"].keys()
        ):
            protein_dict["proteins"][row["uniprot_id"]]["positions"][row["position"]][
                "modifications"
            ].append(
                {
                    "modification_unimod_name": row["modification_unimod_name"],
                    "modification_classification": row["classification"],
                    "modification_unimod_id": row["modification_unimod_id"],
                }
            )

    with open("example_json.json", "w+") as f:
        json.dump(protein_dict, f, indent=3)
    return protein_dict


def parse_user_input(user_file, user_flag):
    try:
        json = read_file(user_file, user_flag)
    except TypeError:
        raise TypeError(
            "Your file could not be parsed. Have you selected the right format?"
        )
    return json


if __name__ == "__main__":
    parse_user_input("example_data/example_data.csv", "csv")
