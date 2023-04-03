from setuptools import setup

setup(
    name="ptmvis",
    packages=["ptmvis"],
    include_package_data=True,
    install_requires=[
        "flask",
        "python-dotenv",
        "Flask-Session",
        "numpy",
        "datetime",
        "numpy",
        "brotli",
        "gunicorn",
    ],
)
