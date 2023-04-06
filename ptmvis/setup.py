from setuptools import setup

setup(
    name="ptmvis",
    packages=["ptmvis"],
    include_package_data=True,
    install_requires=[
        "flask",
        "Flask-Session",
        "python-dotenv",
        "numpy",
        "datetime",
        "brotli",
        "gunicorn",
        "psm_utils"
    ],
)
