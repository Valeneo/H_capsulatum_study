#!/usr/bin/env bash
# Index the WU24 reference once. Every aligner/caller needs its own index or dictionary.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$HOME/miniconda3/envs/hcap/bin:$PATH"
REF="$ROOT/reference/WU24.fa"

# bwa-mem2 index: the FM-index (BWT + suffix array) from lectures 19-20. Built once, searched for every read.
bwa-mem2 index "$REF"
# samtools faidx: byte offsets of each chromosome so tools can jump straight to a position.
samtools faidx "$REF"
# GATK sequence dictionary: chromosome names + lengths + checksum; GATK refuses to run without it.
gatk CreateSequenceDictionary -R "$REF" -O "${REF%.fa}.dict"

# Nuclear chromosomes only. The mitochondrion (CP069113.1) is present in hundreds of copies
# per cell, so its read depth is wildly different and it is not a diploid-style nuclear locus.
# Keeping it would distort depth filters and the SNP matrix.
grep '>' "$REF" | grep -v mitochondrion | sed 's/>//;s/ .*//' > "$ROOT/reference/nuclear_chromosomes.list"
cat "$ROOT/reference/nuclear_chromosomes.list"
