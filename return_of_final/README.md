# H. capsulatum comparative genomics pipeline

183 Illumina paired-end WGS isolates (NCBI SRA) -> variants, assemblies, phylogeny, population structure, comparative analysis.
Every stage has a script in `scripts/` and a lesson in `lessons/` explaining how it works and why each choice was made.

Reference: WU24 (GCA_017310585.1). Environments: `conda env hcap` (bwa-mem2, orthofinder, admixture, clipkit); system tools for the rest.

| # | Objective | Scripts | Tools |
|---|---|---|---|
| 1 | Preprocessing / QC | 01_qc_trim.sh, 01b, 01c | vdb-validate, fasterq-dump, fastp |
| 2 | Alignment + variant calling | 00_prepare_reference.sh, 02_align_call.sh, 03_joint_genotype_filter.sh | bwa-mem2, samtools, GATK4 |
| 3 | Assembly + gene prediction | 04_assemble_busco.sh, 05_predict_genes.sh | SPAdes, BUSCO, Augustus |
| 4 | Phylogenomics, population genetics | 06_popgen.sh, 06b_dapc.R | PLINK2, ADMIXTURE, adegenet (DAPC), IQ-TREE2 |
| 5 | Comparative analysis + metadata integration | | OrthoFinder, R |

Verify any finished stage: `python3 scripts/99_verify.py 1|2|3`.

## Running it step by step

Every step is resumable per sample: stop whenever you need the machine (`bash stop.sh`), start the same command later and it
continues where it left off. Run from this folder. Check each step with `python3 scripts/99_verify.py N` before the next.

| Step | Command | Roughly | Check |
|---|---|---|---|
| 1. QC + trimming (183 samples) | `bash scripts/01c_run_all.sh 3 6` | ~4 h total | `99_verify.py 1` |
| 2. Align + call variants | `HC_XMX=2g bash scripts/02c_run_all.sh 4 4` | 6-10 h | `99_verify.py 2` |
| 3. Joint genotyping + filtering | `bash scripts/03_joint_genotype_filter.sh 8` | 1-2 h | final SNP count printed |
| 4. Population genetics | `bash scripts/06_popgen.sh 8` then `Rscript scripts/06b_dapc.R` | <1 h | |
| 5. Assemblies + BUSCO, full depth, ONE at a time (reads are regenerated from the .sra just in time and deleted after) | `ls results/01_qc \| grep SRR \| xargs -P 1 -I{} bash scripts/04_assemble_busco.sh {} 16` | days | `99_verify.py 3` |
| 6. Gene prediction | `ls results/04_assembly \| grep SRR \| xargs -P 10 -I{} bash scripts/05_predict_genes.sh {}` | 1-2 h | genes ~9,000 per isolate |
| 7. Virulence profiling | `ls results/04_assembly \| grep SRR \| xargs -P 10 -I{} bash scripts/08_virulence.sh {} 2; python3 scripts/08b_virulence_summary.py` | <1 h | |
| 8. Phylogeny | `bash scripts/07_phylogeny.sh 20 1000` | hours | |

`run_pipeline.sh` chains all of it unattended, for when the machine is free for a day or more.

**After a shutdown or crash:** run `bash stop.sh` first (clears half-written samples and large temp files), then restart the step.
Disk: trimmed reads total ~290 GB; once assemblies (step 5) are done they can be deleted (`rm -r results/01_qc/*/*.fastq.gz`).

**Disk policy (space is scarce, time is not):** nothing is capped or subsampled. Trimmed reads are deleted as soon as a sample's gVCF exists and are
regenerated from the `.sra` just before assembly (`REGEN=1`), then deleted again. Set `KEEP_TRIMMED=1` to keep them. Guards wait whenever free
space would fall below `DISK_RESERVE_GB` (default 20) plus the step's working space.

**CPU cap:** `bash cap_cpu.sh 50` pins the running pipeline (and everything it spawns) to the highest-numbered logical CPUs, i.e. 10 of 20
(the 8 E-cores + 1 P-core on this i7-12700H), leaving the fast P-cores free. Re-run it after starting a step from a new shell. `cap_cpu.sh 100` removes the cap.

**Throttling:** `bash throttle.sh keep RUN1 RUN2` keeps only the named step-2 samples running and freezes the rest (nothing lost, frozen memory gets swapped out);
`bash throttle.sh full` unfreezes everything. Combine with `bash cap_cpu.sh 15` (3 CPUs) for a low-impact trickle, or `bash cap_cpu.sh 50` for the normal cap.
