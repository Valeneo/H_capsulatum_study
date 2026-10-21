#!/usr/bin/env bash
# Download and checksum-verify the WU24 reference (NCBI GCA_017310585.1, complete genome, NAm1).
# WU24 is the NAm1 reference of Sepulveda et al. 2017; 7 nuclear chromosomes + mitochondrion, 32.5 Mb.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; mkdir -p "$ROOT/reference"; cd "$ROOT/reference"
B=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/310/585/GCA_017310585.1_ASM1731058v1
F=GCA_017310585.1_ASM1731058v1_genomic
curl -sSfL -O $B/$F.fna.gz -O $B/$F.gff.gz -O $B/md5checksums.txt
grep -E " \./$F\.(fna|gff)\.gz$" md5checksums.txt | md5sum -c -
gunzip -kf $F.fna.gz && mv -f $F.fna WU24.fa
