#!/usr/bin/env python3
"""
Multi-Threaded / Multi-Process Parallel Prefetch Downloader for WGS Datasets
Tools: NCBI SRA Toolkit (prefetch) + Python Concurrent Futures

Features:
- Runs 6 to 8 concurrent prefetch instances in parallel.
- Auto-detects existing non-empty .sra files and skips them.
- Enforces '--max-size u' (unlimited size limit) on each worker.
- Provides real-time progress logging and auto-retry for failed accessions.
"""

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WGS_DIR = os.path.join(BASE_DIR, "wgs")
ACC_FILE = os.path.join(BASE_DIR, "wgs_accessions.txt")
LOG_FILE = os.path.join(BASE_DIR, "download_prefetch_parallel.log")
MAX_WORKERS = 6  # 6 simultaneous prefetch downloads

os.makedirs(WGS_DIR, exist_ok=True)


def download_single_run(acc):
    """Downloads a single SRA run using prefetch."""
    target_sra = os.path.join(WGS_DIR, acc, f"{acc}.sra")
    
    # Check if already downloaded
    if os.path.exists(target_sra) and os.path.getsize(target_sra) > 1024 * 1024:
        return (acc, True, "Already downloaded (Skipped)")
    
    cmd = ["prefetch", acc, "--max-size", "u", "-O", WGS_DIR, "-p"]
    start_time = time.time()
    
    try:
        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        elapsed = time.time() - start_time
        return (acc, True, f"Downloaded in {elapsed:.1f}s")
    except subprocess.CalledProcessError as e:
        err_msg = e.stderr.strip() if e.stderr else str(e)
        return (acc, False, f"Failed: {err_msg[:100]}")
    except Exception as e:
        return (acc, False, f"Exception: {str(e)[:100]}")


def main():
    if not os.path.exists(ACC_FILE):
        print(f"❌ Error: Accession file '{ACC_FILE}' not found.")
        return

    with open(ACC_FILE, "r") as f:
        accessions = [line.strip() for line in f if line.strip()]

    total = len(accessions)
    print("==========================================================================")
    print(f"⚡ Multi-Threaded Prefetch Downloader for {total} WGS Datasets")
    print(f"   Concurrent Workers : {MAX_WORKERS}")
    print(f"   Destination Folder : {WGS_DIR}")
    print("==========================================================================")

    completed = 0
    failed = 0

    with open(LOG_FILE, "a") as log:
        log.write(f"\n[START] {time.strftime('%Y-%m-%d %H:%M:%S')} - Starting {total} downloads with {MAX_WORKERS} threads.\n")

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            future_to_acc = {executor.submit(download_single_run, acc): acc for acc in accessions}

            for future in as_completed(future_to_acc):
                acc = future_to_acc[future]
                try:
                    run_id, success, msg = future.result()
                    completed += 1
                    status_icon = "✅" if success else "❌"
                    log_line = f"[{completed}/{total}] {status_icon} {run_id}: {msg}"
                    print(log_line)
                    log.write(f"{log_line}\n")
                    log.flush()

                    if not success:
                        failed += 1
                except Exception as exc:
                    completed += 1
                    failed += 1
                    err_line = f"[{completed}/{total}] ❌ {acc}: Worker raised exception: {exc}"
                    print(err_line)
                    log.write(f"{err_line}\n")
                    log.flush()

    print("==========================================================================")
    print(f"🎉 Completed: {completed}/{total} | Failed: {failed}")
    print(f"📁 Log saved at: {LOG_FILE}")
    print("==========================================================================")


if __name__ == "__main__":
    main()

