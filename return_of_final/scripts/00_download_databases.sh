#!/usr/bin/env bash
# Databases: BUSCO Onygenales lineage (assembly completeness + phylogenetic marker genes) and Pfam-A (protein domains for HMMER).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"; mkdir -p reference/databases
busco --download onygenales_odb10 --download_path "$ROOT/busco_downloads"
curl -sSfL -C - -o reference/databases/Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
