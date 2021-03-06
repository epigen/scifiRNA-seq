#! /usr/bin/env python

import sys


def parse_requirements(req_file):
    reqs = open(req_file).read().strip().split("\n")
    reqs = [r for r in reqs if not r.startswith("#")]
    return [r for r in reqs if "#egg=" not in r]


# take care of extra required modules depending on Python version
extra = {}
try:
    from setuptools import setup, find_packages

    if sys.version_info < (2, 7):
        extra["install_requires"] = ["argparse"]
    if sys.version_info >= (3,):
        extra["use_2to3"] = True
except ImportError:
    from distutils.core import setup

    if sys.version_info < (2, 7):
        extra["dependencies"] = ["argparse"]

# Requirements
requirements = parse_requirements("requirements.txt")
# requirements_test = parse_requirements(
#     "requirements/requirements.test.txt")
# requirements_docs = parse_requirements(
#     "requirements/requirements.docs.txt")

long_description = open("README.md").read()


# setup
setup(
    name="scifi",
    packages=find_packages(),
    use_scm_version={
        "write_to": "scifi/_version.py",
        "write_to_template": '__version__ = "{version}"\n',
    },
    entry_points={
        "console_scripts": [
            "scifi = scifi.pipeline:main",
            "scifi-summarizer = scifi.scripts.summarizer:main",
        ]
    },
    description="A processing pipeline for scifiRNA-seq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: "
        "GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics, sequencing, ngs, ngs analysis, "
    "single cell, RNA-seq, scRNA-seq",
    url="https://github.com/epigen/scifiRNA-seq",
    project_urls={
        "Bug Tracker": "https://github.com/epigen/scifiRNA-seq/issues",
        "Documentation": "https://scifiRNA-seq.readthedocs.io",
        "Source Code": "https://github.com/epigen/scifiRNA-seq",
    },
    author=u"Andre Rendeiro",
    author_email="andre.rendeiro@pm.me",
    license="GPL3",
    setup_requires=["setuptools_scm"],
    install_requires=requirements,
    # tests_require=requirements_test,
    extras_require={
        # "testing": requirements_test,
        # "docs": requirements_docs
    },
    package_data={"scifi": ["config/default.yaml", "requirements.txt"]},
    **extra
)
