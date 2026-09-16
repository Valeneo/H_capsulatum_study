#!/bin/bash
set -euo pipefail

# Directory configuration
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WGS_DIR="${BASE_DIR}/wgs"
ACC_FILE="${BASE_DIR}/wgs_accessions.txt"
LOG_FILE="${BASE_DIR}/download_wgs.log"

mkdir -p "${WGS_DIR}"

echo "=========================================================="
echo "Histoplasma capsulatum WGS Downloader"
echo "Target directory: ${WGS_DIR}"
echo "Accession list:   ${ACC_FILE}"
echo "Log file:         ${LOG_FILE}"
echo "=========================================================="

total_runs=$(wc -l < "${ACC_FILE}" | tr -d ' ')
current=0
failed=0

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting download of ${total_runs} WGS accessions..." | tee -a "${LOG_FILE}"

while IFS= read -r acc || [ -n "$acc" ]; do
    acc=$(echo "$acc" | tr -d '\r\n[:space:]')
    [ -z "$acc" ] && continue
    
    current=$((current + 1))
    target_sra="${WGS_DIR}/${acc}/${acc}.sra"
    
    # Check if file already exists and is non-empty
    if [ -f "${target_sra}" ] && [ -s "${target_sra}" ]; then
        echo "[${current}/${total_runs}] ${acc} already downloaded. Skipping."
        continue
    fi
    
    echo "----------------------------------------------------------"
    echo "[${current}/${total_runs}] Downloading ${acc}..."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading ${acc}..." >> "${LOG_FILE}"
    
    if prefetch "${acc}" --max-size u -O "${WGS_DIR}" -p >> "${LOG_FILE}" 2>&1; then
        echo "[${current}/${total_runs}] SUCCESS: ${acc} downloaded."
    else
        echo "[${current}/${total_runs}] WARNING: Failed to download ${acc}. See ${LOG_FILE}."
        failed=$((failed + 1))
    fi
done < "${ACC_FILE}"

echo "=========================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished downloads!" | tee -a "${LOG_FILE}"
echo "Total: ${total_runs} | Failed: ${failed}" | tee -a "${LOG_FILE}"
echo "=========================================================="
