from ptmvis import app
from ptmvis.backend import utils
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
    response = []
    if MODIFICATIONS_DATA in session:
        for protein_identifier in session[MODIFICATIONS_DATA]["proteins"]:
            modified_positions = []
            modifications = []
            for modified_position in session[MODIFICATIONS_DATA]["proteins"][
                protein_identifier
            ]["positions"]:
                modified_positions.append(modified_position)
                for modification in session[MODIFICATIONS_DATA]["proteins"][
                    protein_identifier
                ]["positions"][modified_position]["modifications"]:
                    modifications.append(modification["modification_unimod_name"])
            modifications = list(set(modifications))
            protein_entry = {
                "id": protein_identifier,
                "modified_positions": len(modified_positions),
                "unique_modifications": len(modifications),
            }
            response.append(protein_entry)
            print(
                requests.get(
                    "https://rest.uniprot.org/uniprotkb/" + protein_identifier + ".json"
                ).content.decode()
            )
    return response


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
