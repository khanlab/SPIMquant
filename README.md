# SPIMquant

[![Documentation Status](https://readthedocs.org/projects/spimquant/badge/?version=latest)](https://spimquant.readthedocs.io/en/latest/?badge=latest)

SPIMquant is a Snakebids app for quantitative analysis of SPIM (lightsheet) brain data. It performs automated nonlinear template registration and quantification of pathology from SPIM microscopy datasets.

# Table of Contents
 - [Introduction](#introduction)
 - [Installation](#installation)
 - [Usage](#usage)
 - [Contributing](#contributing)
 - [License](#license)

# Introduction
 SPIMquant is designed to process and analyze lightsheet microscopy data. It provides automated tools for nonlinear template registration and quantification of various pathologies.

# Installation

## Hardware Requirements
 - Sufficient memory (at least 16G of memory) for `greedy` diffeomorphic registration.

## Software Requirements
 - A Linux machine with Singularity or Apptainer installed. Alternatively, the following libraries should be available on the command line for a Windows machine:
   - `itk-snap` (includes `c3d`)
   - `greedy` command line tool
   - `python3` with the environment set up according to `pyproject.toml`

## Steps
 1. Clone the repository:
    ```bash
    git clone https://github.com/khanlab/SPIMquant.git
    ```
 2. Install the package:
    ```bash
    pip install -e git+https://github.com/khanlab/spimquant#egg=spimquant
    ```

# Usage
 1. Perform a dry run:
    ```bash
    spimquant /path/to/bids/dir /path/to/output/dir participant -np --use-apptainer
    ```
 2. Run the app using all cores:
    ```bash
    spimquant /path/to/bids/dir /path/to/output/dir participant --cores all --use-apptainer
    ```

# Contributing
 We welcome contributions! Please refer to the [contributing guidelines](CONTRIBUTING.md) for more details on how to contribute.

# License
 This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
