#!/usr/bin/env bash
# Run the whole pipeline end to end. Every step is resumable: finished samples are skipped, so after a crash or
# reboot just run this again. Stops at the first step that fails its checks, and says which one.
# Usage: bash run_pipeline.sh            (logs: logs/pipeline.log)
set -uo pipefail
cd "$(dirname "$0")"
LOG=logs/pipeline.log; mkdir -p logs

say()  { echo "[$(date +%F\ %H:%M)] $*" | tee -a "$LOG"; }
step() { local name="$1"; shift; say "START $name"; if "$@" >>"logs/step_${name}.log" 2>&1; then say "OK    $name"; else say "FAIL  $name  (see logs/step_${name}.log)"; exit 1; fi; }
check() { python3 scripts/99_verify.py "$1" | tail -1 | tee -a "$LOG"; }

# wait for a stage-1 batch that may already be running in another shell
while pgrep -f '[s]cripts/01c_run_all.sh' >/dev/null; do sleep 30; done

# --- Objective 1: QC --------------------------------------------------------------------------------------------
step stage1_qc        bash scripts/01c_run_all.sh 3 6;           check 1
# --- Objective 2: alignment, variants -------------------------------------------------------------------------
[[ -s reference/WU24.fa.bwt.2bit.64 ]] || step prepare_reference bash scripts/00_prepare_reference.sh
step stage2_align     env HC_XMX=2g bash scripts/02c_run_all.sh 4 4;  check 2
step stage2_joint     bash scripts/03_joint_genotype_filter.sh 8
# --- Objective 3: assemblies, genes ---------------------------------------------------------------------------
step stage3_assembly  bash -c 'ls results/01_qc/*/*.fastp.json | xargs -n1 dirname | xargs -n1 basename | xargs -P 2 -I{} bash scripts/04_assemble_busco.sh {} 8';  check 3
step stage3_genes     bash -c 'ls results/04_assembly/*/busco_summary.txt | xargs -n1 dirname | xargs -n1 basename | xargs -P 10 -I{} bash scripts/05_predict_genes.sh {}'
# --- Objective 4: phylogeny, population genetics, virulence -----------------------------------------------------
step popgen           bash scripts/06_popgen.sh 8
step dapc             Rscript scripts/06b_dapc.R
step virulence        bash -c 'ls results/04_assembly/*/busco_summary.txt | xargs -n1 dirname | xargs -n1 basename | xargs -P 10 -I{} bash scripts/08_virulence.sh {} 2 && python3 scripts/08b_virulence_summary.py'
step phylogeny        bash scripts/07_phylogeny.sh 20 1000
say "PIPELINE COMPLETE"
