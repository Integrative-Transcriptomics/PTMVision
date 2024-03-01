from ptmvis import app
from ptmvis.ptmvis_api import api_parameters
from flask import render_template


@app.route("/")
def get_landing_page():
    return render_template("index.html", api_parameters=api_parameters)

@app.route("/about")
def get_about_page():
    return render_template("about.html")
