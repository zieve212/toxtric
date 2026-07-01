"""
setup.py
Packaging for ToxReadout so it installs as a real `toxreadout` command.

Install (from the project root, inside your virtualenv):
    pip install -e .

Then run:
    toxreadout --compound-name "caffeine" --protein-name "adenosine A2A receptor"
"""

from setuptools import setup, find_packages

setup(
    name="toxreadout",
    version="0.1.0",
    description="Compound-protein interaction toxicology readout from "
                "PubChem BioAssay, UniProt/NCBI, and PubMed data.",
    packages=find_packages(exclude=["tests", "tests.*"]),
    python_requires=">=3.10",
    install_requires=[
        "rdkit",
        "biopython",
        "requests",
        "click",
        "rich",
    ],
    extras_require={
        "test": ["pytest"],
    },
    entry_points={
        "console_scripts": [
            "toxreadout=toxreadout.cli:main",
        ],
    },
)
