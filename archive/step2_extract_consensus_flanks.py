#!/usr/bin/env python3
"""
Step 2: Merge, Deduplicate, and Extract 300bp Flanking Sequences for Histoplasma capsulatum
Tools: BioPython (SeqIO)

Biology & Logic:
- MISA detects exact repeats while TRF detects repeats with indels/mismatches.
- We intersect coordinates (overlap >= 15 bp) between MISA and TRF to eliminate redundancy
  and retain high-confidence STR loci.
- For each valid consensus STR locus, we extract exactly 300 bp 5' (upstream) and 300 bp 3' (downstream)
  flanks. Loci near contig ends lacking complete 300 bp flanks are filtered out.
- Normalized chromosome mapping ensures seamless matching across FASTA, MISA, and TRF formats.
"""

import os
import glob
from Bio import SeqIO

GENOME_DIR = "genomes"
MISA_OUT_DIR = "misa/outputs"
TRF_OUT_DIR = "trf/outputs"
OUTPUT_DIR = "flanking_sequences"
FLANK_LENGTH = 300  # 300 bp upstream and downstream

os.makedirs(OUTPUT_DIR, exist_ok=True)


def get_fasta_id_map(fasta_path):
    """Loads a mapping of all FASTA contig IDs and full header descriptions."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    fa_ids = [r.id for r in records]
    return fa_ids, {r.id: r.seq for r in records}


def match_to_fasta_id(raw_chrom, fa_ids):
    """Maps a raw chromosome string from MISA or TRF to the true FASTA contig ID."""
    if raw_chrom in fa_ids:
        return raw_chrom
    # First token match
    first_token = raw_chrom.split()[0]
    if first_token in fa_ids:
        return first_token
    # Prefix match (for MISA underscore replacements)
    for fid in fa_ids:
        if raw_chrom == fid or raw_chrom.startswith(fid + "_") or raw_chrom.startswith(fid):
            return fid
    return None


def parse_misa(misa_path, fa_ids):
    """Parses MISA output file (.misa) and maps to FASTA IDs."""
    loci = []
    if not os.path.exists(misa_path):
        return loci
    with open(misa_path, "r") as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) < 7:
            continue
        raw_chrom, ssr_nr, ssr_type, motif, size, start, end = parts[:7]
        chrom = match_to_fasta_id(raw_chrom, fa_ids)
        if chrom:
            loci.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "motif": motif
            })
    return loci


def parse_trf(trf_path, fa_ids):
    """Parses TRF output file (.dat) and maps to FASTA IDs."""
    loci = []
    if not os.path.exists(trf_path):
        return loci
    current_chrom = None
    with open(trf_path, "r") as f:
        for line in f:
            line_str = line.strip()
            if line_str.startswith("Sequence:"):
                # Extract description after 'Sequence:'
                seq_desc = line_str.replace("Sequence:", "").strip()
                current_chrom = match_to_fasta_id(seq_desc, fa_ids)
                continue
            parts = line_str.split()
            if len(parts) == 15 and current_chrom and parts[0].isdigit():
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    period = int(parts[2])
                    if period in [2, 3, 4, 5, 6]:
                        loci.append({
                            "chrom": current_chrom,
                            "start": start,
                            "end": end,
                            "period": period
                        })
                except ValueError:
                    continue
    return loci


def find_genome_files():
    """Finds all available genome assemblies."""
    datasets = []
    fasta_extensions = ["*.fasta", "*.fna", "*.fa"]
    genome_files = []
    for ext in fasta_extensions:
        genome_files.extend(glob.glob(os.path.join(GENOME_DIR, ext)))

    # Filter out non-FASTA files
    genome_files = [f for f in genome_files if not f.endswith('.misa') and not f.endswith('.statistics') and os.path.getsize(f) > 1024*1024]

    for fasta in sorted(genome_files):
        base_name = os.path.splitext(os.path.basename(fasta))[0]
        misa_file = os.path.join(MISA_OUT_DIR, f"{base_name}.misa")
        trf_file = os.path.join(TRF_OUT_DIR, f"{base_name}.dat")
        datasets.append({
            "name": base_name,
            "fasta": fasta,
            "misa": misa_file,
            "trf": trf_file
        })
    return datasets


def main():
    datasets = find_genome_files()
    if not datasets:
        print(f"❌ No genome assemblies found in '{GENOME_DIR}'.")
        return

    print("==========================================================================")
    print(f"🧬 Histoplasma capsulatum Flank Extraction - Stage 2 (Consensus 300bp)")
    print("==========================================================================\n")

    for dataset in datasets:
        print(f"--- Processing Consensus STRs for: {dataset['name']} ---")

        if not os.path.exists(dataset["fasta"]):
            print(f"  ❌ Missing FASTA: {dataset['fasta']}. Skipping.")
            continue
        if not os.path.exists(dataset["misa"]):
            print(f"  ⚠️ Missing MISA: {dataset['misa']}. Skipping.")
            continue
        if not os.path.exists(dataset["trf"]):
            print(f"  ⚠️ Missing TRF: {dataset['trf']}. Skipping.")
            continue

        fa_ids, genome_dict = get_fasta_id_map(dataset["fasta"])
        misa_loci = parse_misa(dataset["misa"], fa_ids)
        trf_loci = parse_trf(dataset["trf"], fa_ids)

        # Index TRF by contig ID for instant coordinate lookups
        trf_by_chrom = {}
        for t in trf_loci:
            trf_by_chrom.setdefault(t["chrom"], []).append(t)

        consensus_loci = []
        consensus_count = 0

        for m in misa_loci:
            chrom = m["chrom"]
            if chrom not in trf_by_chrom:
                continue

            for t in trf_by_chrom[chrom]:
                overlap_start = max(m["start"], t["start"])
                overlap_end = min(m["end"], t["end"])

                if overlap_start < overlap_end:
                    overlap_len = overlap_end - overlap_start + 1
                    # Overlap >= 15 bp between MISA and TRF
                    if overlap_len >= 15:
                        consensus_count += 1
                        consensus_loci.append({
                            "chrom": chrom,
                            "start": m["start"],
                            "end": m["end"],
                            "motif": m["motif"],
                            "locus_id": f"STR_HC_{dataset['name']}_{consensus_count:05d}"
                        })
                        break

        print(f"  ✅ High-confidence consensus loci (MISA + TRF overlap >= 15bp): {len(consensus_loci):,}")

        # Extract 300 bp flanks
        output_fasta_path = os.path.join(OUTPUT_DIR, f"{dataset['name']}_consensus_flanks.fasta")
        extracted_count = 0

        with open(output_fasta_path, "w") as out_fa:
            for locus in consensus_loci:
                chrom = locus["chrom"]
                start_idx = locus["start"] - 1
                end_idx = locus["end"]

                if chrom in genome_dict:
                    chrom_seq = genome_dict[chrom]
                    seq_len = len(chrom_seq)

                    left_start = max(0, start_idx - FLANK_LENGTH)
                    right_end = min(seq_len, end_idx + FLANK_LENGTH)

                    # Strict check: Require complete 300 bp on BOTH flanks
                    if (start_idx - left_start) == FLANK_LENGTH and (right_end - end_idx) == FLANK_LENGTH:
                        left_flank = chrom_seq[left_start:start_idx]
                        repeat_region = chrom_seq[start_idx:end_idx]
                        right_flank = chrom_seq[end_idx:right_end]

                        full_record = str(left_flank + repeat_region + right_flank)
                        header = f"{locus['locus_id']}|{chrom}|{locus['start']}-{locus['end']}|{locus['motif']}"

                        out_fa.write(f">{header}\n{full_record}\n")
                        extracted_count += 1

        print(f"  📦 Extracted {extracted_count:,} validated 300bp-flanked STR records to:\n     {output_fasta_path}\n")

    print("==========================================================================")
    print("🎉 STAGE 2 FINISHED SUCCESSFULLY!")
    print("==========================================================================")


if __name__ == "__main__":
    main()
