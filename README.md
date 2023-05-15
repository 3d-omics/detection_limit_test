# Snakemake workflow: `Detection limit test`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for assessing detection limit from laser-microdissected samples.


## Usage

### References
Reference genomes (or their soft links) must be placed in the `resources/reference` directory. Ensure that the relative paths of the references in the `config/features.yaml` file match with the actual files.

### Input data
Input data information must be provided in a tabular format through the `config/samples.tsv` file.
