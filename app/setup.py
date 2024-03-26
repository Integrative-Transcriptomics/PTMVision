from setuptools import setup

setup(
    name="ptmvision",
    packages=["ptmvision"],
    version="1.1",
    author="Caroline Jachmann,Simon Hackl",
    author_email="caroline.jachmann@ugent.be,simon.hackl@uni-tuebingen.de",
    description="PTMVision is an interactive web platform for effortless exploratory visual analysis of post translational modifications (PTMs).",
    include_package_data=True,
    install_requires=[
        "flask",
        "Flask-Session",
        "python-dotenv",
        "numpy",
        "datetime",
        "brotli",
        "gunicorn",
        "psm_utils",
        "biopython",
        "unimod-mapper",
        "scipy",
    ],
)