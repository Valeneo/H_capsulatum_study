#!/usr/bin/env bash
# Download + checksum-verify the proteomes used as reference (WU24) and outgroup (Blastomyces dermatitidis ER-3, RefSeq).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; mkdir -p "$ROOT/reference/proteomes"; cd "$ROOT/reference/proteomes"
get() { curl -sSfL -O "$1/$2_protein.faa.gz" && curl -sSfL -o md5.tmp "$1/md5checksums.txt" && grep " \./$2_protein.faa.gz" md5.tmp | md5sum -c - && gunzip -kf "$2_protein.faa.gz" && mv "$2_protein.faa" "$3"; rm -f md5.tmp; }
get https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/017/310/585/GCA_017310585.1_ASM1731058v1 GCA_017310585.1_ASM1731058v1 WU24_ref.faa
get https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/525/GCF_000003525.1_BD_ER3_V1 GCF_000003525.1_BD_ER3_V1 Blastomyces_dermatitidis_outgroup.faa
