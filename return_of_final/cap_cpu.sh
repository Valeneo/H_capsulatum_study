#!/usr/bin/env bash
# Cap the pipeline to a share of the CPU by pinning it (and everything it starts) to the HIGHEST-numbered logical CPUs.
# On this i7-12700H those are the slower E-cores (12-19) first, so the fast P-cores stay free for you.
#   bash cap_cpu.sh 50     pipeline may use 10 of 20 logical CPUs (10-19)
#   bash cap_cpu.sh 100    remove the cap (all CPUs)
# Works on processes that are already running; new processes inherit the pinning from their (pinned) parents.
# Re-run it after starting any new step from a fresh shell.
cd "$(dirname "$0")"
PCT="${1:-50}"; TOTAL=$(nproc); N=$(( (TOTAL * PCT + 99) / 100 )); (( N > TOTAL )) && N=$TOTAL
SET="$((TOTAL - N))-$((TOTAL - 1))"
descendants() { for p in "$@"; do echo "$p"; descendants $(pgrep -P "$p"); done; }
ROOTS=$(pgrep -f '[s]cripts/0[0-9]|[r]un_pipeline|[x]args -P|[0-9]c_run_all')
PIDS=$(descendants $ROOTS | sort -un)
for p in $PIDS; do taskset -a -cp "$SET" "$p" >/dev/null 2>&1; done
echo "cap ${PCT}%: pinned $(echo $PIDS | wc -w) processes to logical CPUs $SET ($N of $TOTAL)"
