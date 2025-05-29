fermo_core
=========

[![DOI](https://zenodo.org/badge/671395100.svg)](https://zenodo.org/doi/10.5281/zenodo.11259126) [![PyPI version](https://badge.fury.io/py/fermo_core.svg)](https://badge.fury.io/py/fermo_core)

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

For general information on FERMO, see the [FERMO Metabolomics GitHub Organization](https://github.com/fermo-metabolomics) page.

## Documentation

The high-level documentation can be found [HERE](https://fermo-metabolomics.github.io/fermo_docs/).

For the Sphinx-generated library documentation, see [HERE](https://fermo-metabolomics.github.io/fermo_core/).

## System Requirements

### Hardware requirements

`fermo_core` can be run on a standard computer and does not have any special requirements.

### Software requirements

#### OS Requirements

Local installation of this package is only supported for Linux (tested on Ubuntu 20.04 and 22.04).

#### Python dependencies

Dependencies including exact versions are specified in the [pyproject.toml](./pyproject.toml) file.

## Installation Guide

Depending on the chosen installation option, install time can vary from a few seconds (`uv`) to a few minutes (`conda`).

### With `pip` from PyPI

*Nota bene: we recommend installing `fermo_core` in some kind of virtual environment.*

```commandline
pip install fermo_core
```

Once installed, run as specified in [Run with `pip`](#run-with-pip)

### With `uv` from GitHub

*Assumes that `uv` is installed*

```commandline
git clone git@github.com:fermo-metabolomics/fermo_core.git
cd fermo_core
uv sync
```

Once installed, run as specified in [Run with `uv`](#run-with-hatch)

### With `conda` from GitHub

*Assumes that `conda` is installed*

```commandline
conda create --name fermo_core python=3.11 -y
conda activate fermo_core
git clone git@github.com:fermo-metabolomics/fermo_core.git
cd fermo_core
pip install -e .
```

Once installed, run as specified in [Run with `conda`](#run-with-conda)

## Quick Start

### Running `fermo_core` on your data

As minimal requirement, `fermo_core` takes LC-MS(/MS) metabolomics data, which it can integrate with a range of optional orthogonal data formats.
Compatible formats are described in the [Documentation](https://fermo-metabolomics.github.io/fermo_docs/home/input_output/).

`fermo_core` requires all parameters to be described in a parameters file. This file must follow specifications outlined in the [JSON Schema](fermo_core/config/schema.json) file. For an example, see [case_study_parameters.json](example_data/case_study_parameters.json).

For a more user-friendly version of FERMO, see [FERMO Online](https://fermo.bioinformatics.nl).

### Run with `pip`

```commandline
fermo_core --parameters <your_parameter_file.json>
```

### Run with `uv`

```commandline
uv run fermo_core --parameters <your_parameter_file.json>
```


### Run with `conda`

```commandline
python3 fermo_core/main.py --parameters <your_parameter_file.json>
```


## Demo

### Overview

To demonstrate the functionality of `fermo_core`, we provide an [example dataset](./example_data) sourced from [this publication](https://doi.org/10.1021/acs.jnatprod.0c00807).
It describes a set of extracts from strains belonging to the bacterial genus *Planomonospora* grown in the same condition, showing differential antibiotic activity against *Staphylococcus aureus*.
`fermo_core` can be used to investigate and prioritize the phenotype-associated and differentially abundant molecular features.
Calculation of the `Phenotype Score` results in the selection of a group of molecular features annotated as siomycins, thiopeptides with known anti-Gram positive antibiotic activity.

Details on the experimental conditions can be found in the [Wiki](https://github.com/fermo-metabolomics/fermo_core/wiki/Demo-example-files-methods).

### Setup

All parameters and input data are specified in a [parameters.json](example_data/case_study_parameters.json) file.

### Run the example

Execution time is hardware-dependent but usually takes only a few minutes. 
On a machine running Ubuntu 22.04 with Intel® Core™ i5-7200U CPU @ 2.50GHz x 4 with 8 GiB Memory, execution time was 104 seconds.

#### Run command

*Nota bene: the exact command depends on the type of installation as specified in the [Installation Guide](#installation-guide).*

```commandline
uv run fermo_core --parameters ./example_data/case_study_parameters.json
```

### Results and Interpretation

After successful completion of the run, all results files can be found in `example_data/results`.
The `out.fermo.session.json` file can be visualized in [FERMO online](https://fermo.bioinformatics.nl/) or inspected with a text viewer.
Alternatively, the spreadsheet format file `out.fermo.abbrev.csv` can be inspected.

Antibiotic activity is attributable to the thiopeptide siomycin and congeners (e.g. feature ID `83`).

## Attribution

### License

`fermo_core` is an open source tool licensed under the MIT license (see [LICENSE](LICENSE.md)).

### Publications

See [CITATION.cff](CITATION.cff) or [FERMO online](https://fermo.bioinformatics.nl/) for information on citing `fermo_core`.


## For Developers

*Nota bene: for details on how to contribute to the FERMO project, please refer to [CONTRIBUTING](CONTRIBUTING.md).*

### Development

Instructions for setting up a development environment.

#### Package Installation

*Assumes that `hatch` is installed*

```commandline
hatch env create dev
hatch run dev:pre-commit install
hatch run dev:pytest --run_slow
```

### Documentation

Instructions on setting up and deploying the automated documentation found [here](http://fermo-metabolomics.github.io/fermo_core/).

The documentation rebuilds and deploys automatically on every release using GitHub Actions.

#### Package Installation

*Assumes that `hatch` is installed*

```commandline
hatch env create doc
hatch run doc:sphinx-build -b html docs/source/ docs/_build
```


