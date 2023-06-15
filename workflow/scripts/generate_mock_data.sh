#!/usr/bin/env bash

set -euo pipefail

# Slice
samtools faidx \
    resources/reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz \
    22:20000000-20010000 \
| cut -f 1 -d ":" \
| bgzip -l 9 > resources/reference/human_22_sub.fa.gz

samtools faidx \
    resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz \
    39:100000-110000 \
| cut -f 1 -d ":" \
| bgzip -l 9 > resources/reference/chicken_39_sub.fa.gz
samtools faidx \
    resources/reference/mags.fa.gz \
    NZ_AP023069.1:1-10000 \
    NC_015714.1:1-10000 \
| cut -f 1 -d ":" \
| bgzip -l 9 > resources/reference/mags_sub.fa.gz


# Simulate
mkdir --parents resources/reads/

wgsim \
    -N 1000 \
    -S 1 \
    resources/reference/human_22_sub.fa.gz \
    >(gzip -9 > resources/reads/human1_1.fq.gz) \
    >(gzip -9 > resources/reads/human1_2.fq.gz) \
> resources/reads/human1.log 2>&1

wgsim \
    -N 1000 \
    -S 2 \
    resources/reference/human_22_sub.fa.gz \
    >(gzip -9 > resources/reads/human2_1.fq.gz) \
    >(gzip -9 > resources/reads/human2_2.fq.gz) \
> resources/reads/human2.log 2>&1

wgsim \
    -N 1000 \
    -S 1 \
    resources/reference/chicken_39_sub.fa.gz \
    >(gzip -9 > resources/reads/chicken1_1.fq.gz) \
    >(gzip -9 > resources/reads/chicken1_2.fq.gz) \
> resources/reads/chicken1.log 2>&1

wgsim \
    -N 1000 \
    -S 2 \
    resources/reference/chicken_39_sub.fa.gz \
    >(gzip -9 > resources/reads/chicken2_1.fq.gz) \
    >(gzip -9 > resources/reads/chicken2_2.fq.gz) \
> resources/reads/chicken2.log 2>&1


samtools faidx resources/reference/mags_sub.fa.gz

wgsim \
    -N 1000 \
    -S 1 \
    resources/reference/mags_sub.fa.gz \
    >(gzip -9 > resources/reads/mags1_1.fq.gz) \
    >(gzip -9 > resources/reads/mags1_2.fq.gz) \
> resources/reads/mags1.log 2>&1

wgsim \
    -N 1000 \
    -S 2 \
    resources/reference/mags_sub.fa.gz \
    >(gzip -9 > resources/reads/mags2_1.fq.gz) \
    >(gzip -9 > resources/reads/mags2_2.fq.gz) \
> resources/reads/mags2.log 2>&1


# Merge
gzip -dc \
    resources/reads/human1_1.fq.gz \
    resources/reads/chicken1_1.fq.gz \
    resources/reads/mags1_1.fq.gz \
| bgzip -l 9 > resources/reads/sample1_1.fq.gz

gzip -dc \
    resources/reads/human1_2.fq.gz \
    resources/reads/chicken1_2.fq.gz \
    resources/reads/mags1_2.fq.gz \
| bgzip -l 9 > resources/reads/sample1_2.fq.gz

gzip -dc \
    resources/reads/human2_1.fq.gz \
    resources/reads/chicken2_1.fq.gz \
    resources/reads/mags2_1.fq.gz \
| bgzip -l 9 > resources/reads/sample2_1.fq.gz

gzip -dc \
    resources/reads/human2_2.fq.gz \
    resources/reads/chicken2_2.fq.gz \
    resources/reads/mags2_2.fq.gz \
| bgzip -l 9 > resources/reads/sample2_2.fq.gz
