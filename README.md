# Snakemake workflow: `detection limit test`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/3d-omics/bioinfo_detection_limit_test/workflows/Tests/badge.svg?branch=devel)](https://github.com/3d-omics/bioinfo_detection_limit_test/actions?query=branch%devel+workflow%3ATests)


A Snakemake workflow for assessing detection limit from laser-microdissected samples.

## Usage

0. Requirements
   1.  [`git-lfs`](https://git-lfs.com/)
   2.  [`miniconda`](https://docs.conda.io/en/latest/miniconda.html) / [`mamba`](https://mamba.readthedocs.io)
   3.  [`snakemake`](snakemake.readthedocs.io/)

1. Clone the repository
Clone the repository, and set it as the working directory.

```
git clone --recursive https://github.com/3d-omics/detection_limit_test.git
cd detection_limit_test
```

2. Run the pipeline with the test data (takes 5 minutes)
```
snakemake \
    --use-conda \
    --conda-frontend mamba \
    -j 8
```

3. Edit the following files:
   1. `config/samples.tsv`: the control file with the sequencing libraries and their location.
   2. `config/features.tsv`: the references against which to map the libraries: human, chicken / pig, MAG catalogue.
   3. `config/params.tsv`: parameters for every program. The defaults are reasonable.

4. Run the pipeline on a cluster with slurm and go for a walk:

```
./run_slurm
```

## Brief description
1. Trim reads and remove adaptors with `fastp`
2. Map to human, chicken / pig, mag catalogue:
   1. Map to the reference
   2. Extract the reads that have one of both ends unmapped
   3. Map those unmapped reads to the next reference
3. Generate MAG statistics with
   1. CoverM
   2. SingleM
   3. Nonpareil
4. Generate lots of reports in the `reports/` folder
