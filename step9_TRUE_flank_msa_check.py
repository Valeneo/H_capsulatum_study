import os
import subprocess
import re
from Bio import SeqIO

META_FILE = "results/step2_conserved_loci/conserved_loci_metadata.txt"
POLY_FILE = "results/step4_polymorphism/polymorphic_loci_ranked.tsv"
MSA_DIR = "results/stage5_msa_validation"
os.makedirs(MSA_DIR, exist_ok=True)

TRUE_VALIDATED_FILE = os.path.join(MSA_DIR, "true_msa_validated_loci.tsv")

# Databases paths for fetching real sequence blocks
DB_PATHS = {
    "b68": "blast/databases/barcode68_db",
    "cuti": "blast/databases/cuti_ref_db",
    "cybja": "blast/databases/cybja_ref_db"
}

def fetch_sequence_from_db(db_path, seq_id, start, end):
    """Uses blastdbcmd to pull real sequence fragments from a local BLAST DB."""
    # Handle flipped orientations from BLAST coordinates
    s_min, s_max = min(int(start), int(end)), max(int(start), int(end))
    
    # Run blastdbcmd
    cmd = f"blastdbcmd -db {db_path} -entry '{seq_id}' -range {s_min}-{s_max} -outfmt %s"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    sequence = res.stdout.strip().upper()
    
    # If aligned on reverse strand, return reverse complement
    if int(start) > int(end):
        from Bio.Seq import Seq
        sequence = str(Seq(sequence).reverse_complement())
    return sequence

# Load top 50 highly ranked polymorphic IDs to test
target_loci = []
with open(POLY_FILE, "r") as f:
    lines = f.readlines()
    for line in lines[1:51]:  # Let's test the top 50 first to be safe
        target_loci.append(line.split("\t")[0])

# Parse coordinates from metadata mapping file
metadata = {}
with open(META_FILE, "r") as f:
    for line in f.readlines()[1:]:
        p = line.strip().split("\t")
        if p:
            metadata[p[0]] = p

print(f"🔬 Starting TRUE MSA validation on the top {len(target_loci)} polymorphic candidates...")

validated_count = 0

with open(TRUE_VALIDATED_FILE, "w") as out:
    out.write("Locus_ID\tTrue_MSA_Identity\tStatus\n")
    
    for locus_id in target_loci:
        if locus_id not in metadata:
            continue
        
        p = metadata[locus_id]
        
        # 1. Gather Left Flank sequence for Barcode 67 from its header context
        # Our original pipeline extracted 300bp left flank
        # We can simulate/fetch the exact left flank block
        header_parts = locus_id.split("|")
        
        # 2. Fetch the true corresponding regions across B68, Cuti, and Cybja
        # We look at the target start coordinate and grab the 300bp upstream flanking anchor window
        try:
            # Fetch real sequence slices from your actual assembly files
            seq_b68  = fetch_sequence_from_db(DB_PATHS["b68"], p[4], int(p[5])-300, p[5])
            seq_cuti = fetch_sequence_from_db(DB_PATHS["cuti"], p[7], int(p[8])-300, p[8])
            seq_cybja= fetch_sequence_from_db(DB_PATHS["cybja"], p[10], int(p[11])-300, p[11])
            
            # Fallback length check
            if not seq_b68 or not seq_cuti or not seq_cybja:
                continue
                
            temp_fasta = os.path.join(MSA_DIR, "real_temp.fasta")
            temp_aligned = os.path.join(MSA_DIR, "real_aligned.fasta")
            
            with open(temp_fasta, "w") as tmp:
                tmp.write(f">b68\n{seq_b68}\n")
                tmp.write(f">cuti\n{seq_cuti}\n")
                tmp.write(f">cybja\n{seq_cybja}\n")
                
            # Run MAFFT over the actual independent sequences
            subprocess.run(f"mafft --quiet {temp_fasta} > {temp_aligned}", shell=True)
            
            # Calculate alignment score
            records = list(SeqIO.parse(temp_aligned, "fasta"))
            seq_len = len(records[0].seq)
            matches = 0
            for col in range(seq_len):
                chars = [r.seq[col].upper() for r in records if r.seq[col] != '-']
                if len(chars) == 3 and len(set(chars)) == 1:
                    matches += 1
            
            identity = (matches / seq_len) * 100
            
            status = "PASSED" if identity >= 90.0 else "DISCARDED"
            out.write(f"{locus_id}\t{identity:.2f}%\t{status}\n")
            if identity >= 90.0:
                validated_count += 1
                
        except Exception as e:
            # Handle out-of-bounds near edge of contigs gracefully
            continue

print(f"\n✨ Real MSA Check Finished!")
print(f"  - Verified {validated_count} loci with genuine biological conservation >= 90% across strains.")
print(f"📊 True verification log written to: {TRUE_VALIDATED_FILE}")
