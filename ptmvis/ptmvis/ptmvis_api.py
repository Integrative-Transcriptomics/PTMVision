from ptmvis import app
from flask import request, session, render_template, send_file
from flask_session import Session
from datetime import timedelta
from dotenv import load_dotenv
import json, zlib, os, subprocess, shutil, numpy, brotli, random, string

""" Load variables from local file system. """
load_dotenv()

""" Set session configuration parameters. """
app.config["SESSION_TYPE"] = "filesystem"
app.config["SESSION_KEY_PREFIX"] = "ptmvis:"
app.config["SESSION_COOKIE_NAME"] = "ptmvis_session"
app.config["SESSION_FILE_DIR"] = "./ptmvis/session/"
app.config["SESSION_USE_SIGNER"] = False
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_FILE_THRESHOLD"] = 15000
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
