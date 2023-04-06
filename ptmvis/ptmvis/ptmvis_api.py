from ptmvis import app
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
app.config["SESSION_COOKIE_NAME"] = "ptmvis_session"
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


@app.route("/process_search_engine_output", methods=["POST"])
def process_search_engine_output():
    """
    TODO
    """
    response = "Ok"
    try:
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        content_type = json_request_data["contentType"]
        df = pd.read_csv(StringIO(json_request_data["content"]))
        print(df)
    except Exception as e:
        print(e)
    finally:
        return response
