from ptmvision import app
from ptmvision.api import DEBUG
from flask import render_template


@app.route("/")
def get_landing_page():
    return render_template("index.html")

@app.route("/about")
def get_about_page():
    return render_template("about.html")
