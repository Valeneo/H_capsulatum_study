#!/usr/bin/env bash
# Safely stop whatever pipeline step is running, then clean half-written outputs so the step can be resumed later.
# Finished samples are never touched. Usage: bash stop.sh
cd "$(dirname "$0")"
PAT='[s]cripts/0[0-9]|[r]un_pipeline|[x]args -P|[f]astp -i|[f]asterq-dump|[s]pades|[b]wa-mem2|[G]enomicsDBImport|[H]aplotypeCaller|[a]ugustus|[i]qtree2|[m]afft|[b]usco'
PIDS=$(pgrep -f "$PAT"); [ -n "$PIDS" ] && kill -9 $PIDS 2>/dev/null; sleep 2
rm -rf results/01_qc/.tmp_* results/02_variants/*/.sorttmp results/02_variants/*/*.bam results/02_variants/*/*.bam.bai
# a sample is only finished if its LAST-written file exists; anything else is partial
for d in results/01_qc/SRR*/;       do r=$(basename $d); [ -s $d$r.fastp.json ]       || rm -rf "$d"; done
for d in results/02_variants/SRR*/; do r=$(basename $d); [ -s $d$r.g.vcf.gz.tbi ]     || rm -rf "$d"; done
for d in results/04_assembly/SRR*/; do r=$(basename $d); [ -s $d/busco_summary.txt ] || rm -rf "$d"; done
echo "stopped. complete: stage1=$(ls results/01_qc/*/*.fastp.json 2>/dev/null | wc -l) stage2=$(ls results/02_variants/*/*.g.vcf.gz.tbi 2>/dev/null | wc -l) stage3=$(ls results/04_assembly/*/busco_summary.txt 2>/dev/null | wc -l)"
