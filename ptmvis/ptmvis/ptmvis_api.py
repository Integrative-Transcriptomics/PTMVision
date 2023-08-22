from ptmvis import app
from ptmvis.backend import utils
from ptmvis.backend.uniprot_id_mapping import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)
from flask import request, session, render_template, send_file
from flask_session import Session
from dotenv import load_dotenv
from io import StringIO
import json, zlib, os, requests

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
        # Try to map UniProt identifiers to gene names.
        gene_name_mapping_response = _map_uniprot_identifiers(
            protein_identifiers, "Gene_Name"
        )
        gene_name_mapping = {}
        # Collect mapping results into dictionary.
        for entry in gene_name_mapping_response["results"]:
            gene_name_mapping[entry["from"]] = entry["to"]
        if "failedIds" in gene_name_mapping_response:
            for failed_identifier in gene_name_mapping_response["failedIds"]:
                gene_name_mapping[failed_identifier] = "N/A"
        # Construct entry per protein in input data.
        for protein_identifier in protein_identifiers:
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
                "name": gene_name_mapping[protein_identifier],
                "modified_positions": len(modified_positions),
                "unique_modifications": len(modifications),
                "modifications": "$".join(modifications),
            }
            protein_entries.append(protein_entry)
    return protein_entries


@app.route("/get_protein_data", methods=["POST"])
def get_protein_data():
    """
    TODO
    """
    response = {"status": "Ok", "ptms": None, "contacts": None}
    try:
        json_request_data = request_to_json(request.data)
        structure = None
        if json_request_data["opt_pdb_text"] != None:
            structure = utils.parse_structure(json_request_data["opt_pdb_text"])
        else:
            structure = utils.get_structure(json_request_data["uniprot_id"])
        if structure != None:
            df = session[MODIFICATIONS_DATA]
            response["ptms"] = df[
                df["uniprot_id"] == json_request_data["uniprot_id"]
            ].to_csv()
            response["contacts"] = utils.get_contacts(
                utils.get_distance_matrix(structure),
                int(json_request_data["distance_cutoff"]),
            )
            response["sequence"] = utils.get_sequence_from_structure(structure)
    except Exception as e:
        response["status"] = "Failed: " + str(e)
    finally:
        return response


def _map_uniprot_identifiers(identifiers, target_db):
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db=target_db, ids=identifiers
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    return results
