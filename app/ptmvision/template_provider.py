from ptmvision import app
from ptmvision.api import DEBUG
from flask import render_template


@app.route("/")
def route_index():
    return render_template("index.html")

@app.route("/usage")
def route_usage_guide():
    return render_template("usage_guide.html")

@app.route("/about")
def route_about():
    return render_template("about.html")

@app.route("/legal")
def route_legal_information():
    return render_template("legal_information.html")

@app.route("/faq")
def route_faq():
    return render_template("faq.html")
