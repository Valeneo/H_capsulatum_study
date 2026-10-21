# Lesson 1: Preprocessing for quality control

Objective: *To perform preprocessing for quality control.*
Script: `scripts/01_qc_trim.sh` (one sample), `scripts/01b_qc_summary.py` (all samples).
Worked example below is real output from SRR25758805.

## The logic: why this step exists

Everything downstream (alignment, variant calling, assembly) assumes the reads are
trustworthy. Raw reads contain three kinds of junk, all from your Lecture 10 notes:

| Problem | Where it comes from | What it breaks |
|---|---|---|
| Adapter sequence in the read | Insert shorter than the read, so the read runs through into the adapter | The aligner cannot place the read (adapter has no genomic match) |
| Low-quality bases (3' end decays) | Chemistry degrades along the read | False SNPs |
| Very short / mostly-bad reads after trimming | Consequence of the two above | Ambiguous, multi-mapping reads |

Preprocessing removes these *before* any analysis, because errors that enter early
get amplified later and cannot be seen afterwards.

## Step by step, and why each step is there

**1. `vdb-validate`: is the .sra file intact?**
Your `wgs/` folder has leftover `.aria2`/`.lock` files from interrupted downloads.
A truncated .sra can still be opened and give *fewer reads with no error*. That would
silently bias one sample. Validate first, fail loudly.

**2. `fasterq-dump --split-files`: .sra to FASTQ**
.sra is NCBI's compressed archive; no aligner reads it. `--split-files` makes two files
per sample, `_1` and `_2`: the two ends of each DNA fragment (paired-end).
The FASTQ is ~4-5x larger than the .sra, so we delete it as soon as fastp is done
(`trap ... rm -rf "$TMP"`).

**3. `fastp`: trim and QC in one pass**

| Option | Meaning | Reason |
|---|---|---|
| `--detect_adapter_for_pe` | Finds adapters by looking where mate 1 and mate 2 overlap | No adapter file to get wrong across 183 isolates from different labs |
| `--cut_right -W 4 -M 20` | Slide a 4-base window from the 3' end, cut when mean quality drops under 20 | Removes the decayed tail, keeps the good 5' part of the read (the "quality trimming" of Lecture 10) |
| `-q 20 -u 30` | A base below Q20 is "bad"; drop the read if over 30% of bases are bad | Removes reads that are bad throughout |
| `-l 50` | Drop reads under 50 bp after trimming | Short reads map to many places |
| `--trim_poly_g` | Remove runs of G at the read end | On two-colour chemistry (NovaSeq/NextSeq), "no signal" is read as G |
| no `--dedup` | Deliberately not removing duplicates here | Duplicates are found more accurately later from mapping coordinates (`samtools markdup`) |

**Why Q20 and not Q30?** The detailed flowchart said "Q >= 30". As a trimming cutoff
that is too aggressive: it throws away usable bases and lowers coverage, which hurts
variant calling more than it helps. Q20 (1 error in 100) for trimming, then judge the
result by how much of the data is Q30 afterwards.

## Reading the result (SRR25758805)

| Metric | Before | After | Meaning |
|---|---|---|---|
| Reads | 2,012,486 | 1,908,940 | 94.9% kept |
| Bases | 480 Mb | 434 Mb | ~9.6% of bases removed |
| Q30 | 93.8% | 98.3% | The point of the step: cleaner data |
| GC | 45.9% | 45.5% | Barely moves, so trimming did not bias composition |
| Duplication | | 0.4% | Low, good library complexity |
| Reads dropped as too short | | 87,674 | Largest loss, mostly tails trimmed to under 50 bp |
| Adapter-trimmed reads | | 3,660 | Few: inserts were long enough (peak insert size 399 bp) |

Coverage estimate: 434 Mb / ~35 Mb genome is about 12x. That is on the low side and is
worth remembering when we set variant-calling filters. (Genome size to be confirmed
from the reference we pick.)

## Things to know before the next step

- **Not done yet: contamination screening (Kraken2).** There is no Kraken2 database on
  this machine; it needs a download (the standard DB is ~50-70 GB) and the disk is 94% full.
- **Raw/clean data does not fit on disk.** The 183 clean FASTQ.gz will be ~250 GB; 90 GB is free.
