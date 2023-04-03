from flask import Flask

app = Flask(__name__)

from ptmvis import ptmvis_api
from ptmvis import template_provider

if __name__ == "__main__":
    app.run()
