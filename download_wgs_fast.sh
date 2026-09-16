#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WGS_DIR="${BASE_DIR}/wgs"
ACC_FILE="${BASE_DIR}/wgs_accessions.txt"
LOG_FILE="${BASE_DIR}/download_fast.log"
INPUT_FILE="${BASE_DIR}/aria2_input.txt"

mkdir -p "${WGS_DIR}"

echo "=========================================================================="
echo "⚡ High-Performance Multi-Threaded WGS Downloader (via aria2c & AWS Mirror)"
echo "Target Directory   : ${WGS_DIR}"
echo "Accession List     : ${ACC_FILE}"
echo "Log File           : ${LOG_FILE}"
echo "=========================================================================="

total_runs=$(wc -l < "${ACC_FILE}" | tr -d ' ')
> "${INPUT_FILE}"

while IFS= read -r acc || [ -n "$acc" ]; do
    acc=$(echo "$acc" | tr -d '\r\n[:space:]')
    [ -z "$acc" ] && continue
    
    target_dir="${WGS_DIR}/${acc}"
    target_file="${target_dir}/${acc}.sra"
    
    # Check if already downloaded and non-empty
    if [ -f "${target_file}" ] && [ -s "${target_file}" ]; then
        continue
    fi
    
    mkdir -p "${target_dir}"
    # AWS SRA Open Data Direct Mirror + NCBI Fallback Mirror
    echo "https://sra-pub-run-odp.s3.amazonaws.com/sra/${acc}/${acc}" >> "${INPUT_FILE}"
    echo "  dir=${target_dir}" >> "${INPUT_FILE}"
    echo "  out=${acc}.sra" >> "${INPUT_FILE}"
done < "${ACC_FILE}"

remaining=$(grep -c "^https" "${INPUT_FILE}" || true)
echo "Total Runs in Queue: ${remaining} / ${total_runs}"

if [ "${remaining}" -eq 0 ]; then
    echo "🎉 All 183 WGS datasets are already downloaded in ${WGS_DIR}!"
    exit 0
fi

echo "🚀 Launching 8 concurrent downloads with 16 parallel connections each..."
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting aria2c parallel download..." >> "${LOG_FILE}"

# Optimized for 20-core systems: 8 concurrent files, 16 connections per file
aria2c \
    --input-file="${INPUT_FILE}" \
    --max-concurrent-downloads=8 \
    --split=16 \
    --min-split-size=5M \
    --max-connection-per-server=16 \
    --continue=true \
    --auto-file-renaming=false \
    --allow-overwrite=true \
    --summary-interval=5 \
    --log="${LOG_FILE}" \
    --log-level=notice
