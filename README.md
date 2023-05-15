# Snakemake workflow: `Detection limit test`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for assessing detection limit from laser-microdissected samples.

## Usage

### Clone repository
Clone the repository, and set it as the working directory.

```
git clone https://github.com/3d-omics/detection_limit_test.git
cd detection_limit_test
```

### References
Reference genomes (or their soft links) must be placed in the `resources/reference` directory. Ensure that the relative paths of the references in the `config/features.yaml` file match with the actual files.

```
mkdir resources
mkdir resources/reference
cd resources/reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
cd ../..
```

### Input data
Input data information must be provided in a tabular format through the `config/samples.tsv` file.

```
mkdir resources/reads
cd resources/reads/
wget https://sid.erda.dk/share_redirect/EG29XsTN5s/D056cJa1_EKDL220017891-1A_HKJMJDSX5_L2_1.fq.gz
wget https://sid.erda.dk/share_redirect/EG29XsTN5s/D056cJa1_EKDL220017891-1A_HKJMJDSX5_L2_2.fq.gz
wget https://sid.erda.dk/share_redirect/EG29XsTN5s/D018aJa1_EKDL220017891-1A_HKJMJDSX5_L2_1.fq.gz
wget https://sid.erda.dk/share_redirect/EG29XsTN5s/D018aJa1_EKDL220017891-1A_HKJMJDSX5_L2_2.fq.gz
cd ../..
```

### Run snakemake

```
module load snakemake/7.20.0 mamba/1.3.1
snakemake --dry-run --use-conda --jobs 8 all
```
