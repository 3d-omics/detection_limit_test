#!/usr/bin/env bash
set -euo pipefail

# rename all fastas in the folder to fasta_filename^header

for f in *.fa; do
    sed -i "s/^>/>${f}\^/ ; s/.fa//" "$f";
done
