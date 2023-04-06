import re
import pandas as pd
from Bio import SeqIO
from io import BytesIO
from unimod_mapper import UnimodMapper

mapper = UnimodMapper()

def get_ptm_position_in_protein(mod, protein_start):
    position_ptm_in_peptide = int(re.search(".+?(?=\D)", mod)[0]) 
    position_ptm_in_protein = protein_start + position_ptm_in_peptide - 1
    return position_ptm_in_protein

def map_mass_to_unimod(mods):
    modstr = ""
    for mod in mods.split(","):
        pos_aa = mod.split("(")[0]

        mass = re.search(r'\(.*?\)', mod)[0]
        mass = float(mass.replace("(", "").replace(")", ""))
        mapped = mapper.mass_to_names(mass, decimals=4) # for unimod ID: mapper.mass_to_ids(mass, decimals = 4)

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
    protein_df = pd.DataFrame(zipped, columns=["uniprot_id","position", "modification_unimod_name"])
    protein_df["modification_unimod_id"] = protein_df["modification_unimod_name"].apply(lambda x: id_from_name(x))

    return protein_df

def read_msfragger(file):
    # read file and pick relevant rows and columns
    print("Reading MSFragger file...")
    pept_mods = pd.read_csv(file, sep = "\t")
    pept_mods = pept_mods[["Modified Peptide", "Peptide", "Protein ID", "Assigned Modifications", "Protein Start"]]
    pept_mods = pept_mods.dropna(subset=["Modified Peptide"])
    pept_mods = pept_mods.drop_duplicates()

    # try to assign modifications to mass shifts
    print("Assigning modifications to mass shifts...")
    pept_mods["Assigned Modifications parsed"] = pept_mods["Assigned Modifications"].apply(lambda x: map_mass_to_unimod(x)) 
    
    # rewrite to protein level modifications
    print("Mapping modification sites onto proteins...")
    prot_mods = from_psm_to_protein(pept_mods)

    return prot_mods

def read_ionbot(file):
    df = pd.read_csv(file)
    df = df[["uniprot_id", "unexpected_modification", "position"]]
    df["modification_unimod_name"] = df["unexpected_modification"].apply(lambda x: x.split("]")[1].split("[")[0].lower())
    df["modification_unimod_id"] = df["unexpected_modification"].apply(lambda x: x.split("]")[0][1:])
    df.drop(columns=["unexpected_modification"], inplace=True)
    return df

def read_mod_csv(file):
    df = pd.read_csv(file)
    return df


def read_file(file, flag):
    if flag == "ionbot":
        df = read_ionbot(file)
    elif flag == "msfragger":
        df = read_msfragger(file)
    elif flag == "csv":
        df = read_mod_csv(file)
    return df


def parse_user_input(user_file, user_flag):
    try:
        df = read_file(user_file, user_flag)
    except TypeError:
        print("Your file could not be parsed. Have you selected the right format?")
    return df