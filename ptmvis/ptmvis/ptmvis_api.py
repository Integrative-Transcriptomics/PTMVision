from ptmvis import app
from ptmvis.backend import utils
from ptmvis.backend.uniprot_id_mapping import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)
from flask import request, session, send_file
from flask_session import Session
from dotenv import load_dotenv
from io import StringIO
from itertools import combinations
from math import ceil
from copy import deepcopy
import json, zlib, os

""" Load variables from local file system. """
load_dotenv()

""" Set session configuration parameters. """
app.config["SESSION_TYPE"] = "filesystem"
app.config["SESSION_FILE_DIR"] = "./ptmvis/session/"
app.config["SESSION_KEY_PREFIX"] = "ptmvis:"
app.config["SESSION_USE_SIGNER"] = False
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_FILE_THRESHOLD"] = 1000
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # Limit content lengths to 200 MB.

""" Set API parameters. """
api_parameters = {"URL": os.getenv("URL"), "DEBUG": os.getenv("DEBUG")}

""" Start session """
Session(app)

""" Definition of session keys """
MODIFICATIONS_DATA = "modifications_data_frame"
AVAILABLE_PROTEINS = "available_proteins"


def request_to_json(request_content: any) -> object:
    """Decompress zlib encoded JSON request into Python object.

    Uses UTF-8 decoding.

    Parameters
    ----------
    request_content : ReadableBuffer
        zlib encoded JSON request content.

    Returns
    -------
    object
        Python object yielding the content of the JSON request.
    """
    inflated_request_data = zlib.decompress(request_content)
    json_string_request_data = inflated_request_data.decode("utf8")
    return json.loads(json_string_request_data)


@app.route("/example_data", methods=["GET"])
def example_data():
    """Route to download example CSV format data."""
    return send_file("./static/resources/example_data.csv", as_attachment=True)


@app.route("/process_search_engine_output", methods=["POST"])
def process_search_engine_output():
    """Route to process CSV data or SAGE, ionbot or MSFragger search engine output."""
    response = "Ok"
    try:
        json_request_data = request_to_json(request.data)
        content_type = json_request_data["contentType"]
        json_user_data = utils.parse_user_input(
            StringIO(json_request_data["content"]), content_type
        )
        session[MODIFICATIONS_DATA] = json_user_data
    except TypeError as e:
        response = "Failed: " + str(e)
    finally:
        return response


@app.route("/get_available_proteins", methods=["GET"])
def get_available_proteins():
    """Route to retrieve all available proteins of the session as a JSON."""
    protein_entries = []
    if MODIFICATIONS_DATA in session:
        protein_identifiers = [_ for _ in session[MODIFICATIONS_DATA]["proteins"]]
        with open(
            "./ptmvis/static/resources/uniprot_id_mapping.json", "r"
        ) as protein_name_mapping_in:
            protein_name_mapping = json.load(protein_name_mapping_in)
        # Construct entry per protein in input data.
        for protein_identifier in protein_identifiers:
            if protein_identifier in protein_name_mapping:
                protein_name = protein_name_mapping[protein_identifier]
            else:
                protein_name = "N/A"
            position_modification_data = session[MODIFICATIONS_DATA]["proteins"][
                protein_identifier
            ]["positions"]
            modified_positions = []
            modifications = []
            for modified_position in position_modification_data:
                modified_positions.append(int(modified_position))
                for modification in position_modification_data[modified_position][
                    "modifications"
                ]:
                    modifications.append(modification["modification_unimod_name"])
            modified_positions = sorted(modified_positions)
            modifications = list(set(modifications))
            protein_entry = {
                "id": protein_identifier,
                "name": protein_name,
                "modified_positions": len(modified_positions),
                "unique_modifications": len(modifications),
                "modifications": "$".join(modifications),
            }
            protein_entries.append(protein_entry)
    return protein_entries


@app.route("/get_modifications_graph", methods=["GET"])
def get_modifications_graph():
    """Route to retrieve an ECharts graph definition representing all modifications of a dataset."""
    nodes = {}
    links = {}
    modification_occurrence = {}
    if MODIFICATIONS_DATA in session:
        # Construct full modifications graph; Nodes represent modificiations, Links represent common occurrence.
        protein_identifiers = [_ for _ in session[MODIFICATIONS_DATA]["proteins"]]
        for protein_identifier in protein_identifiers:
            position_modification_data = session[MODIFICATIONS_DATA]["proteins"][
                protein_identifier
            ]["positions"]
            for position in position_modification_data:
                per_position_modifications = []
                for modification in position_modification_data[position][
                    "modifications"
                ]:
                    modification_name = modification["modification_unimod_name"]
                    per_position_modifications.append(modification_name)
                    if not modification_name in nodes:
                        nodes[modification_name] = {
                            "name": modification_name,
                            "value": 1,
                        }
                    else:
                        nodes[modification_name]["value"] += 1
                    if not modification_name in modification_occurrence:
                        modification_occurrence[modification_name] = {
                            protein_identifier
                        }
                    else:
                        modification_occurrence[modification_name].add(
                            protein_identifier
                        )
                for pair in combinations(set(per_position_modifications), 2):
                    try:
                        link_key = frozenset([pair[0], pair[1]])
                        if not link_key in links:
                            links[link_key] = {
                                "source": list(link_key)[0],
                                "target": list(link_key)[1],
                                "value": 1,
                            }
                        else:
                            links[link_key]["value"] += 1
                    except IndexError as e:
                        print(pair)
        # Filter modifications graph for top L and adjust layout settings.
        # (i) Adjust nodes.
        L = 50
        nodes = dict( sorted( nodes.items(), key = lambda n: n[ 1 ][ "value" ], reverse = True )[ :L ] )
        nodes_values = [node["value"] for node in nodes.values()]
        nodes_values_min = min(nodes_values)
        nodes_values_max = max(nodes_values)
        for k, v in nodes.items():
            v["symbolSize"] = max(
                4,
                20 if (nodes_values_max - nodes_values_min) == 0 else
                ceil(
                    (
                        (v["value"] - nodes_values_min)
                        / (nodes_values_max - nodes_values_min)
                    )
                    * 20
                ),
            )
            v["count"] = v["value"]
            v["value"] = round(
                (len(modification_occurrence[k]) / len(protein_identifiers)) * 100
            )
        # (ii) Adjust links.
        links = [
            l for l in links.values() if (l["source"] in nodes and l["target"] in nodes)
        ]
        if len( links ) > 0 :
            links_values = sorted([link["value"] for link in links])
            links_values_min = min(links_values)
            links_values_max = max(links_values)
            for link in links:
                link["lineStyle"] = {
                    "width": max(
                        0.05,
                        2 if (links_values_max - links_values_min) == 0 else
                        (
                            (link["value"] - links_values_min)
                            / (links_values_max - links_values_min)
                        )
                        * 2,
                    )
                }
        # Construct EChart option.
        option = {
            "backgroundColor": "#fbfbfb",
            "animation": False,
            "tooltip": {"position": [5, 45]},
            "visualMap": {
                "bottom": "bottom",
                "left": "center",
                "min": 0,
                "max": 100,
                "color": ["#F26430", "#4350A5"],
                "orient": "horizontal",
                "precision": 1,
                "text": ["", "Frequency in Proteins"],
            },
            "series": [
                {
                    "name": "Modifications Graph",
                    "emphasis": {
                        "focus": "adjacency",
                        "label": {"show": False},
                        "itemStyle": {"shadowColor": "#dc5754", "shadowBlur": 10},
                    },
                    "label": {
                        "show": True,
                        "color": "#333333",
                        "fontWeight": "lighter",
                        "fontSize": 11,
                        "backgroundColor": "rgba(240, 245, 245, 0.6)",
                        "borderRadius": 4

                    },
                    "edgeLabel": {"show": False},
                    "type": "graph",
                    "layout": "circular",
                    "circular": {
                        "rotateLabel": True
                    },
                    "data": list(nodes.values()),
                    "links": links,
                    "roam": True,
                    "lineStyle": {"color": "#333333", "curveness": 0.2},
                    "itemStyle": {},
                }
            ],
        }
    return option


@app.route("/get_protein_data", methods=["POST"])
def get_protein_data():
    response = {"status": "ok", "content": None}
    json_request_data = request_to_json(request.data)
    # Try to fetch PDB format structure for UniProt identifier from AlphaFold database.
    if json_request_data["pdb_text"] != None:
        structure, pdb_text = utils.parse_structure(json_request_data["pdb_text"])
    else:
        if "structure" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] :
            pdb_text = utils._brotli_decompress( session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["structure"] )
            structure = utils.parse_structure( pdb_text )
        else :
            structure, pdb_text = utils.get_structure(json_request_data["uniprot_id"])
            # Extract protein sequence from structure and store it in session data.
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["sequence"] = utils.get_sequence_from_structure(structure)
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "structure" ] = utils._brotly_compress(pdb_text)
    if structure != None:
        # Extract annotations for protein from UniProt.
        if not "annotation" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] :
            annotation = _map_uniprot_identifiers(
                [json_request_data["uniprot_id"]], "UniProtKB"
            )
            # TODO: Clean annotation.
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["annotation"] = { }
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["annotation"] = annotation[ "results" ][ 0 ][ "to" ]

        # Compute contacts from structure and store them in session data.
        if not "contacts" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] :
            session[MODIFICATIONS_DATA]["meta_data"]["distance_cutoff"] = int(json_request_data["cutoff"])
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "contacts" ] = { }
            distance_matrix = utils.get_distance_matrix(structure)
            contacts = utils.get_contacts(
                distance_matrix,
                float(json_request_data["cutoff"]),
            )
            for source_position, contacts_list in contacts.items( ) :
                session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "contacts" ][ source_position ] = [ ]
                for contact_position in contacts_list :
                    session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "contacts" ][ source_position ].append( ( contact_position, distance_matrix[ source_position, contact_position ] ) )
        # Construct response.
        response[ "content" ] = deepcopy( session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] )
        response[ "content" ][ "structure" ] = utils._brotli_decompress( response[ "content" ][ "structure" ] )
    else :
        response = {"status": "Failed: No protein structure available.", "content": None}
    with open('./ptmvis/session/dump.json', 'w+') as f:
        json.dump( session[MODIFICATIONS_DATA], f, indent = 3 )
    return response

def _map_uniprot_identifiers(identifiers, target_db):
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db=target_db, ids=identifiers
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    return results
