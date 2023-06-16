#!/usr/bin/env bash

set -euo pipefail

mkdir --parents resources/reference/

wget \
    --continue \
    --directory-prefix resources/reference/ \
    https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz \
    https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/030/075/GCF_013030075.1_ASM1303007v1/GCF_013030075.1_ASM1303007v1_genomic.fna.gz \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/217/655/GCF_000217655.1_ASM21765v1/GCF_000217655.1_ASM21765v1_genomic.fna.gz

gzip -dc \
    resources/reference/GCF_013030075.1_ASM1303007v1_genomic.fna.gz \
    resources/reference/GCF_000217655.1_ASM21765v1_genomic.fna.gz \
> resources/reference/mags.fa
gzip -d resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz
gzip -d resources/reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

bgzip resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa
bgzip resources/reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa
bgzip resources/reference/mags.fa

samtools faidx resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz
samtools faidx resources/reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
samtools faidx resources/reference/mags.fa.gz
