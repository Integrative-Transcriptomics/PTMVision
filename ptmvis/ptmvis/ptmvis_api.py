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
PROT_STRUC = "protein_structure"


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
    if AVAIL_PROT in session:
        return session[AVAIL_PROT]
    else:
        return []


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


@app.route("/load_protein_structure_into_session", methods=["POST"])
def load_protein_structure_into_session():
    """
    TODO
    """
    response = "Ok"
    try:
        json_request_data = request_to_json(request.data)
        print(json_request_data)
        if json_request_data["opt_pdb_text"] != None:
            session[PROT_STRUC] = utils.parse_structure(
                json_request_data["opt_pdb_text"]
            )
        else:
            session[PROT_STRUC] = utils.get_structure(json_request_data["uniprot_id"])
    except Exception as e:
        response = "Failed: " + str(e)
    finally:
        return response
