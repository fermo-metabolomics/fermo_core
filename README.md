fermo_core
=========

[![DOI](https://zenodo.org/badge/671395100.svg)](https://zenodo.org/doi/10.5281/zenodo.11259126) [![PyPI version](https://badge.fury.io/py/fermo_core.svg)](https://badge.fury.io/py/fermo_core) [![Go Coverage](https://github.com/fermo-metabolomics/fermo_core/wiki/coverage.svg)](https://raw.githack.com/wiki/fermo-metabolomics/fermo_core/coverage.html)

Contents
-----------------
- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Start](#quick-start)
- [Demo](#demo)
- [Attribution](#attribution)
- [For Developers](#for-developers)

## Overview

`fermo_core` is a tool to perform hypothesis-driven prioritization on metabolomics data. 
Besides its use as backend processing tool of the [FERMO dashboard](https://github.com/fermo-metabolomics/fermo), `fermo_core` can be used as a command line interface (CLI) for large-scale data processing and analysis.

This README specifies the use of `fermo_core` as CLI. For a more user-friendly version of FERMO, see [FERMO Online](https://fermo.bioinformatics.nl).

For more information, see the [FERMO Metabolomics](https://github.com/fermo-metabolomics) GitHub Organization page.

## Documentation

## System Requirements

## Installation Guide

### With `pip` from PyPI
- Install `python 3.11.x`
- Create a virtual environment (e.g. venv, conda) and activate it
- Run `pip install fermo_core`
- Once installed, run as specified in [Run with `pip`](#run-with-pip)

### With `hatch` from GitHub
- Install `python 3.11.x`
- Install hatch (e.g. with `pipx install hatch`)
- Download or clone the repository
- (Change into the fermo_core base directory if not already present)
- Run `hatch -v env create`
- Once installed, run as specified in [Run with `hatch`](#run-with-hatch)

### With `conda` from GitHub
- Install conda (e.g. miniconda)
- Create a conda environment with `conda create --name fermo_core python=3.11`
- Activate the conda environment with `conda activate fermo_core`
- Download or clone the repository
- (Change into the fermo_core base directory if not already present)
- Run `pip install -e .`
- Once installed, run as specified in [Run with `conda`](#run-with-conda)

## Quick Start

### Run with `pip`
- `fermo_core --parameters <your_parameter_file.json>`

### Run with `hatch`:
- `hatch run fermo_core --parameters <your_parameter_file.json>`

### Run with `conda`:
- `python3 fermo_core/main.py --parameters <your_parameter_file.json>`

## Demo


`fermo_core` can be used both as a command line interface as well as a library.

All parameters and input data are specified in a `parameters.json` file.
This file must be formatted according to the JSON Schema specified in `fermo_core/config/schema.json`. 
See the example in `example_data/case_study_parameters.json`.

For more information on input and output files, their format, and their purpose, consult the [Documentation](https://fermo-metabolomics.github.io/fermo_docs/home/input_output/).

## Attribution

### License

`fermo_core` is an open source tool licensed under the MIT license (see [LICENSE](LICENSE.md)).

### Publications

See [FERMO online](https://fermo.bioinformatics.nl/) for information on citing `fermo_core`.


## For Developers