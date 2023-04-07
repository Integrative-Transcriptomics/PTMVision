from ptmvis import app
from ptmvis.backend import utils
from flask import request, session, render_template, send_file
from flask_session import Session
from datetime import timedelta
from dotenv import load_dotenv
from io import StringIO
import json, zlib, os, subprocess, shutil, numpy, brotli, random, string
import pandas as pd

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
api_parameters = {
    "EXAMPLE_PARAMETER": "-1",  # Example parameter for debugging.
    "URL": os.getenv("URL"),
}

""" Set DEBUG to true in order to display log as console output. """
DEBUG = bool(int(os.getenv("DEBUG")))

""" Start session """
Session(app)

""" Definition of session keys """
MOD_DF = "modifications_data_frame"
AVAIL_PROT = "available_proteins"


def request_to_json(zlib_req):
    inflated_request_data = zlib.decompress(zlib_req)
    json_string_request_data = inflated_request_data.decode("utf8")
    return json.loads(json_string_request_data)


@app.route("/process_search_engine_output", methods=["POST"])
def process_search_engine_output():
    """
    TODO
    """
    response = "Ok"
    try:
        json_request_data = request_to_json(request.data)
        content_type = json_request_data["contentType"]
        df = utils.parse_user_input(
            StringIO(json_request_data["content"]), content_type
        )
        session[MOD_DF] = df
        session[AVAIL_PROT] = df["uniprot_id"].unique().tolist()
    except TypeError as e:
        response = "Failed: " + str(e)
    finally:
        return response


@app.route("/get_available_proteins", methods=["GET"])
def get_available_proteins():
    """
    TODO
    """
    response = {"prot_names": []}
    if AVAIL_PROT in session:
        response["prot_names"] = session[AVAIL_PROT]
    return response


@app.route("/get_dashboard_data", methods=["GET"])
def get_dashboard_data():
    """
    TODO
    """
    response = "Ok"
    try:
        json_request_data = request_to_json(request.data)
    except Exception as e:
        response = "Failed: " + str(e)
    finally:
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
            df = session[MOD_DF]
            response["ptms"] = df[
                df["uniprot_id"] == json_request_data["uniprot_id"]
            ].to_csv()
            response["contacts"] = utils.get_contacts(
                utils.get_distance_matrix(structure),
                int(json_request_data["distance_cutoff"]),
            )
    except Exception as e:
        response["status"] = "Failed: " + str(e)
    finally:
        return response
