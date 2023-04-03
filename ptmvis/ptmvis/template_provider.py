from ptmvis import app
from ptmvis.ptmvis_api import api_parameters
from flask import render_template


@app.route("/")
def get_template_home():
    return render_template("home.html", api_parameters=api_parameters)
