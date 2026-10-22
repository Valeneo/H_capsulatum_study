#!/usr/bin/env bash
# Throttle step 2 without losing any work.
#   bash throttle.sh 2      keep the 2 samples that have been running longest, FREEZE the rest (progress kept, nothing lost)
#   bash throttle.sh keep RUN1 RUN2   keep the named samples running, freeze the others
#   bash throttle.sh full   unfreeze everything
#   bash throttle.sh status
# Frozen samples keep their batch slot, so no new samples start beyond the ones allowed to run. Frozen memory gets swapped out
# by the kernel as soon as something else needs the RAM.
cd "$(dirname "$0")"
descendants() { for p in "$@"; do echo "$p"; descendants $(pgrep -P "$p"); done; }
samples() { ps -eo pid,etimes,args | awk '/[0]2_align_call.sh SRR/ && !/bash -c/ {print $2, $1, $(NF-1)}' | sort -rn; }
case "${1:-status}" in
  keep)   shift; keep=" $* "   # keep=run names to leave running, freeze every other sample
          samples | while read -r et pid run; do
            if [[ "$keep" == *" $run "* ]]; then kill -CONT $(descendants $pid) 2>/dev/null; echo "run     $run"
            else kill -STOP $(descendants $pid) 2>/dev/null; echo "freeze  $run"; fi
          done ;;
  full)   for pid in $(samples | awk '{print $2}'); do kill -CONT $(descendants $pid) 2>/dev/null; done; echo "all samples running" ;;
  status) samples | while read -r et pid run; do st=$(ps -o stat= -p "$(pgrep -P $pid | head -1)" 2>/dev/null | cut -c1); echo "$run  in for $((et/60)) min  $([[ $st == T ]] && echo FROZEN || echo active)"; done ;;
  [0-9]*) n=$1; i=0
          samples | while read -r et pid run; do i=$((i+1))
            if (( i <= n )); then kill -CONT $(descendants $pid) 2>/dev/null; echo "run     $run ($((et/60)) min in)"
            else kill -STOP $(descendants $pid) 2>/dev/null; echo "freeze  $run ($((et/60)) min in)"; fi
          done ;;
esac
