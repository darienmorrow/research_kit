#! /usr/bin/env python


# --- import -------------------------------------------------------------------------------------


import os

from setuptools import setup, find_packages


# --- define -------------------------------------------------------------------------------------


here = os.path.abspath(os.path.dirname(__file__))


extra_files = []
extra_files.append(os.path.join(here, "CONTRIBUTORS"))
extra_files.append(os.path.join(here, "LICENSE"))
extra_files.append(os.path.join(here, "README.md"))
extra_files.append(os.path.join(here, "research_kit", "VERSION"))


# --- setup --------------------------------------------------------------------------------------


with open(os.path.join(here, "requirements.txt")) as f:
    required = f.read().splitlines()


with open(os.path.join(here, "research_kit", "VERSION")) as version_file:
    version = version_file.read().strip()


setup(
    name="research_kit",
    version=version,
    packages=find_packages(),
    package_data={"": extra_files},
    install_requires=required,
    author="Darien Morrow",
    author_email="darienmorrow@gmail.com",
    license="MIT",
    url="https://github.com/darienmorrow/research_kit",
    keywords="photophysics spectroscopy science",
    entry_points={
        "console_scripts": [
            "dir_PL_work=research_kit.__main__:read_plot_save",
            "dir_hl3=research_kit.__main__:dir_hl3",
        ]
    },
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
    ],
)
