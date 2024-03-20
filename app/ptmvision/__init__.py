from flask import Flask

app = Flask(__name__)

from ptmvision import api
from ptmvision import template_provider

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=5001)
