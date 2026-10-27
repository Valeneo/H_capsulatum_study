# Handoff brief: H. capsulatum WGS pipeline, status and how to take over step 2 on another PC

Written 2026-10-27 on the machine that has run the pipeline so far ("PC-A"). The receiving session ("PC-B") has no link to PC-A,
so everything it needs is in this file and in the repo. Read `PIPELINE_DESIGN.md` (one folder up) first: it holds the design, rules and
shared conventions. Do not change the design or the open decisions listed below; carry out the task in section 3.

## 1. Pipeline status (183 Illumina paired-end isolates, reference WU24)
| Step | Status |
|---|---|
| 1. QC and trimming | DONE for all 183 (fastp). 9 samples kept 48-69% of reads (older runs); they are usable (>=8x after trimming). |
| 2. Alignment + GATK haploid gVCFs | IN PROGRESS: 7 of 183 finished (list in section 2). Slow: 3.5-9.6 h per sample for divergent isolates, because GATK time follows divergence from the reference, not depth. |
| 3. Joint genotyping and filtering | Script tested end to end on 3 isolates (79,804 raw -> 32,413 filtered SNPs). Not run on the full set; needs all gVCFs first. |
| 4. PCA, ADMIXTURE, DAPC | Scripts written, never run (need the step 3 SNP matrix). |
| 5. SPAdes assemblies + BUSCO (full output) | Script tested on 1 isolate (3,578 contigs, 29.7 Mb, BUSCO 93.1% complete). Batch not started; estimated 6-10 days of compute at full depth. |
| 6. Gene prediction (Augustus) | Tested on 1 isolate (8,308 proteins; reference has 9,080). |
| 7. Virulence + mating type (tblastn) | Tested on 1 isolate (9/10 genes found, MAT1-1). |
| 8. Phylogeny (BUSCO genes, MAFFT, ClipKIT, IQ-TREE2) | Scripts written, never run. |

## 2. Step 2 detail
Finished (gVCF + index present, results in `handoff/step2_finished_2026-10-27.tgz`, 140 MB):
SRR13508696 (99.99% mapped, pass), SRR25580010 (91.96%, pass), SRR25580028 (99.84%, pass), SRR25580030 (99.82%, pass),
SRR25758805 (99.87%, pass), SRR8084715 (93.41%, pass), SRR13433424 (86.92%, FAILED the 90% mapping gate: strain G184AR, Panama lineage, divergent from the
NAm1 reference; the gate flags it for review, the gVCF itself is valid; keep it and leave it listed in `logs/needs_review.txt`).
Running on PC-A right now: SRR29905417. PC-A works the queue from the FRONT (`config/step2_queue.txt`, shallow samples first, two slow divergent samples last).
Expected low mappers still queued: SRR20306115 (55.7%, China, GC 38%, possible contamination) and SRR13433402 (H88, African lineage, 80.4%).

## 3. Your task (PC-B): run step 2 on the queue from the BACK, so the two machines meet in the middle
Rules that matter: nothing is depth-capped or subset; a sample is finished only when `results/02_variants/<R>/<R>.g.vcf.gz.tbi` exists (written last);
the script verifies the sample and only then deletes the trimmed reads; keep >=20 GB free disk; edit scripts only by atomic replace; do not use pkill -f with a pattern in your own command.
Hardware note: PC-A is an i7-12700H, 15 GB RAM, one sample at a time; size your parallelism to PC-B (HC_XMX=2g per GATK job, about 3-5 GB RAM per sample at peak, 3-4 CPU threads each).

Setup:
1. `git clone https://github.com/Valeneo/H_capsulatum_study.git` and work in `H_capsulatum_study/return_of_final`.
2. Conda env named `hcap`: `mamba create -n hcap -c conda-forge -c bioconda bwa-mem2 samtools bcftools gatk4 fastp sra-tools seqkit` and make sure `java` (17+) is present. The scripts prepend `~/miniconda3/envs/hcap/bin` to PATH.
3. `bash scripts/00_download_reference.sh` then `bash scripts/00_prepare_reference.sh` (indexes the WU24 reference; about a minute).
4. Unpack the finished results so they are not repeated: `tar xzf step2_finished_2026-10-27.tgz` inside `return_of_final` (copy the file over by USB or cloud; it is not in git).
5. Raw data: each sample is downloaded from SRA into `../wgs/<R>/<R>.sra` (the folder `wgs/` sits next to `return_of_final`).

Run (back of the queue first; skips finished samples; step 2 rebuilds the trimmed reads from the .sra by itself):
```
cd return_of_final
for R in $(tac config/step2_queue.txt); do
  [ -s results/02_variants/$R/$R.g.vcf.gz.tbi ] && continue
  [ -s ../wgs/$R/$R.sra ] || prefetch --max-size 100G -O ../wgs "$R"
  HC_XMX=2g bash scripts/02_align_call.sh "$R" 4 >> logs/stage2_$R.log 2>&1
done
```
Several of these loops can run side by side on different samples (the per-sample lock prevents two jobs on one sample).
If the loop reaches a sample that PC-A already finished, that is fine (it is skipped only if you imported its result); duplicates waste time but are harmless.

## 3b. If PC-B already holds raw data (.sra files) elsewhere
Nothing large needs to move between machines: only the 140 MB results package above. Raw `.sra` files are public (NCBI SRA) and can be reused in place.
- Put or link each file at `../wgs/<R>/<R>.sra` (a symlink is fine): `mkdir -p ../wgs/<R> && ln -s /existing/path/<R>.sra ../wgs/<R>/<R>.sra`. The loop in section 3 skips the download when the file is already there.
- Every file is checked with `vdb-validate` at the start of step 1 and rejected loudly if corrupt. Do not skip this: on PC-A, 7 downloaded files turned out truncated or corrupt
  (SRR8084710, SRR8084708, SRR8347492, SRR8426586, SRR6243645, SRR6243646, SRR6243643) and were replaced by fresh `prefetch` downloads that validated.
- Trimmed reads are temporary (up to 8 GB per sample, rebuilt from the .sra on demand and deleted after verification); only the .sra files and the small results are long-lived.

## 4. Returning results to PC-A
Every few samples or when stopping, pack the new results and move the file to PC-A:
`tar czf step2_pcB_$(date +%F).tgz results/02_variants/*/{*.g.vcf.gz,*.g.vcf.gz.tbi,*.flagstat.txt,*.coverage.tsv,*.verify.log} logs/needs_review.txt`
Verify before relying on a sample: `python3 scripts/99_verify.py 2 --run <R>` (mapping >=90%, depth >=5x, sample name, 7 nuclear chromosomes, haploid genotypes).
A sample that fails the gate keeps its trimmed reads and is listed in `logs/needs_review.txt`; report the list back instead of deleting or rerunning.

## 5. Not your decision (the project owner decides)
- Whether to replace GATK with bcftools in haploid mode for the remaining samples (GATK is far slower than planned; a comparison on finished samples is the proposed first step).
- Whether to add a second reference for divergent lineages (H88 for African, a Panama strain) to reduce mapping bias.
- Deleting `wgs_corrupt/` and the original `.sra` archives.
Report problems and the final list of finished samples; do not start steps 3-8.
