#!/usr/bin/env python3
"""
Step 1: Synchronized Microsatellite (STR) Discovery for Histoplasma capsulatum
Tools: MISA (misa.pl) & Tandem Repeats Finder (TRF)

Biology & Logic:
- MISA detects perfect microsatellites (Di >= 6, Tri >= 5, Tetra >= 5, etc.)
- TRF detects both perfect and imperfect tandem repeats (MinScore=24, MaxPeriod=4)
- This script dynamically discovers all genome assemblies (*.fasta, *.fna, *.fa)
  in the 'genomes' directory and processes them with synchronized thresholds.
- Ultra-verbose output with per-genome telemetry, discovery statistics, and time metrics.
"""

import os
import subprocess
import glob
import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from Bio import SeqIO

GENOME_DIR = "genomes"
MISA_OUT_DIR = "misa/outputs"
MISA_CONFIG = "misa.ini"
TRF_OUT_DIR = "trf/outputs"
MAX_CORES = 15  # Strict cap at 15 CPU cores max

# Ensure output directories exist
os.makedirs(MISA_OUT_DIR, exist_ok=True)
os.makedirs(TRF_OUT_DIR, exist_ok=True)

# Create standard MISA configuration if missing
if not os.path.exists(MISA_CONFIG):
    with open(MISA_CONFIG, "w") as f:
        f.write(
            "definition(min_repeats) 1-10 2-6 3-5 4-5 5-5 6-5\n"
            "interruptions(max_difference_for_2_SSRs) 100\n"
        )
    print(f"⚙️ Created standard MISA configuration: {MISA_CONFIG}")


def count_misa_records(misa_file):
    """Counts detected SSRs in MISA output."""
    if not os.path.exists(misa_file):
        return 0, {}
    ssr_count = 0
    type_counts = {}
    with open(misa_file, "r") as f:
        for line in f.readlines()[1:]:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                ssr_count += 1
                ssr_type = parts[2]  # e.g., p2, p3, p4, c
                type_counts[ssr_type] = type_counts.get(ssr_type, 0) + 1
    return ssr_count, type_counts


def count_trf_records(dat_file):
    """Counts detected tandem repeats in TRF output."""
    if not os.path.exists(dat_file):
        return 0
    tr_count = 0
    with open(dat_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 15 and parts[0].isdigit():
                tr_count += 1
    return tr_count


def process_single_genome(fasta_path):
    """Runs MISA and TRF on a single genome assembly with rich verbose logging."""
    start_time = time.time()
    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    full_name = os.path.basename(fasta_path)
    file_size_mb = os.path.getsize(fasta_path) / (1024 * 1024)

    # Count contigs and total base length
    try:
        contig_count = 0
        total_bp = 0
        for rec in SeqIO.parse(fasta_path, "fasta"):
            contig_count += 1
            total_bp += len(rec.seq)
    except Exception:
        contig_count = "N/A"
        total_bp = int(file_size_mb * 1024 * 1024)

    print(f"\n==========================================================================")
    print(f"🚀 [START] Processing: {full_name}")
    print(f"   📊 Assembly Info : {file_size_mb:.2f} MB | {total_bp:,} bp | {contig_count} contig(s)")
    print(f"==========================================================================")

    # 1. RUN MISA
    print(f"   ▶️ [1/2] Running MISA (Di >= 6, Tri >= 5, Tetra >= 5)...")
    misa_start = time.time()
    misa_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "misa.pl")
    misa_cmd = f"perl {misa_script} {fasta_path} {MISA_CONFIG}" if os.path.exists(misa_script) else f"misa.pl {fasta_path} {MISA_CONFIG}"
    
    misa_dest = os.path.join(MISA_OUT_DIR, f"{base_name}.misa")
    stats_dest = os.path.join(MISA_OUT_DIR, f"{base_name}.statistics")

    try:
        subprocess.run(misa_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        generated_misa = fasta_path + ".misa"
        generated_stats = fasta_path + ".statistics"

        if os.path.exists(generated_misa):
            os.replace(generated_misa, misa_dest)
        if os.path.exists(generated_stats):
            os.replace(generated_stats, stats_dest)

        misa_time = time.time() - misa_start
        total_ssrs, ssr_types = count_misa_records(misa_dest)
        print(f"   ✅ [MISA SUCCESS] {base_name} processed in {misa_time:.2f}s")
        print(f"      - Total Perfect SSRs Found: {total_ssrs:,}")
        if ssr_types:
            type_str = ", ".join(f"{k}:{v}" for k, v in sorted(ssr_types.items()))
            print(f"      - SSR Breakdown: {type_str}")
    except subprocess.CalledProcessError:
        print(f"   ❌ [MISA ERROR] Failed on {base_name}. Ensure misa.pl and perl are available.")
    except Exception as e:
        print(f"   ⚠️ [MISA EXCEPTION] {base_name}: {e}")

    # 2. RUN TRF
    print(f"\n   ▶️ [2/2] Running Tandem Repeats Finder (TRF, Period <= 4, MinScore = 24)...")
    trf_start = time.time()
    trf_cmd = f"trf {fasta_path} 2 7 7 80 10 24 4 -d -h"
    trf_dest = os.path.join(TRF_OUT_DIR, f"{base_name}.dat")

    try:
        subprocess.run(trf_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # TRF writes output to current working directory
        base_dat_pattern = f"{os.path.basename(fasta_path)}*.dat"
        name_dat_pattern = f"{base_name}*.dat"
        dat_files = glob.glob(base_dat_pattern) + glob.glob(name_dat_pattern) + glob.glob(f"{fasta_path}*.dat")
        for dat_file in set(dat_files):
            os.replace(dat_file, trf_dest)

        trf_time = time.time() - trf_start
        total_tr = count_trf_records(trf_dest)
        print(f"   ✅ [TRF SUCCESS] {base_name} processed in {trf_time:.2f}s")
        print(f"      - Total Tandem Repeats Identified: {total_tr:,}")
    except subprocess.CalledProcessError:
        print(f"   ❌ [TRF ERROR] Failed on {base_name}. Ensure trf binary is available.")
    except Exception as e:
        print(f"   ⚠️ [TRF EXCEPTION] {base_name}: {e}")

    elapsed_total = time.time() - start_time
    print(f"--------------------------------------------------------------------------")
    print(f"✨ [COMPLETED] {full_name} finished in {elapsed_total:.2f}s total.")
    print(f"--------------------------------------------------------------------------\n")
    return base_name


def main():
    print("==========================================================================")
    print("🧬 Histoplasma capsulatum STR Discovery Pipeline - Stage 1 (MISA + TRF)")
    print("==========================================================================")

    fasta_extensions = ["*.fasta", "*.fna", "*.fa"]
    genome_files = []
    for ext in fasta_extensions:
        genome_files.extend(glob.glob(os.path.join(GENOME_DIR, ext)))

    # Filter out any non-FASTA or statistics files
    genome_files = [f for f in genome_files if not f.endswith('.misa') and not f.endswith('.statistics') and os.path.getsize(f) > 1024*1024]

    if len(genome_files) == 0:
        print(f"❌ No valid genome FASTA files (>1MB) found in '{GENOME_DIR}'.")
        print(f"   Please verify your assembly files in '{GENOME_DIR}/'.")
        return

    # Enforce strict 15-core cap
    available_cpus = multiprocessing.cpu_count()
    num_cpus = min(MAX_CORES, available_cpus, len(genome_files))
    
    print(f"📂 Genome Directory       : {GENOME_DIR}")
    print(f"📦 Discovered Assemblies  : {len(genome_files)} FASTA file(s)")
    for idx, gf in enumerate(sorted(genome_files), 1):
        sz = os.path.getsize(gf) / (1024 * 1024)
        print(f"   [{idx}] {os.path.basename(gf)} ({sz:.2f} MB)")
    print(f"⚡ CPU Threads Allocation : {num_cpus} threads (Strict cap: {MAX_CORES} cores, {available_cpus - num_cpus} cores kept free)")
    print(f"📁 Output Destinations    : '{MISA_OUT_DIR}/' and '{TRF_OUT_DIR}/'")
    print("==========================================================================\n")

    global_start = time.time()
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        results = list(executor.map(process_single_genome, sorted(genome_files)))

    global_elapsed = time.time() - global_start
    print(f"==========================================================================")
    print(f"🎉 STAGE 1 FINISHED SUCCESSFULLY!")
    print(f"   - Processed Genomes : {len(results)} assemblies")
    print(f"   - Total Clock Time  : {global_elapsed:.2f} seconds ({global_elapsed/60:.2f} minutes)")
    print(f"   - Results Ready For : Step 2 (step2_extract_consensus_flanks.py)")
    print(f"==========================================================================")


if __name__ == "__main__":
    main()
