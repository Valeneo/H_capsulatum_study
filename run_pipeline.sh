#!/usr/bin/env bash
set -euo pipefail

# Make sure scripts are executable
chmod +x run_blast_all.sh step17_external_validation_generator.py validate_primers_cli.py 2>/dev/null || true

python3 run_pipeline.py "$@"

