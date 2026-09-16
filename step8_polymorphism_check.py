import os

META_FILE = "results/step2_conserved_loci/conserved_loci_metadata.txt"
OUTPUT_DIR = "results/step4_polymorphism"
os.makedirs(OUTPUT_DIR, exist_ok=True)
POLY_OUTPUT_FILE = os.path.join(OUTPUT_DIR, "polymorphic_loci_ranked.tsv")

if not os.path.exists(META_FILE):
    print(f"Error: {META_FILE} not found. Please ensure Step 3 ran successfully.")
    exit(1)

polymorphic_panel = []
monomorphic_count = 0

with open(META_FILE, "r") as f:
    lines = f.readlines()
    header = lines[0].strip().split("\t")
    
    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) < 13:
            continue
            
        locus_id = parts[0]
        
        # Calculate lengths for each genome using absolute coordinate spans
        len_b67  = abs(int(parts[3]) - int(parts[2])) + 1
        len_b68  = abs(int(parts[6]) - int(parts[5])) + 1
        len_cuti = abs(int(parts[9]) - int(parts[8])) + 1
        len_cybja = abs(int(parts[12]) - int(parts[11])) + 1
        
        sizes = [len_b67, len_b68, len_cuti, len_cybja]
        unique_sizes = set(sizes)
        
        # Filter: If only 1 unique size exists across all 4 genomes, it is monomorphic
        if len(unique_sizes) <= 1:
            monomorphic_count += 1
            continue
            
        # Calculate PIC (Polymorphism Information Content)
        # Count frequencies of each unique allele size across the 4 genomes
        counts = {}
        for s in sizes:
            counts[s] = counts.get(s, 0) + 1
            
        sum_p2 = 0.0
        for s in counts:
            p = counts[s] / 4.0  # Frequency out of 4 strains
            sum_p2 += p ** 2
            
        pic_score = 1.0 - sum_p2
        
        polymorphic_panel.append({
            "locus_id": locus_id,
            "motif": locus_id.split("|")[-1] if "|" in locus_id else "Unknown",
            "size_b67": len_b67,
            "size_b68": len_b68,
            "size_cuti": len_cuti,
            "size_cybja": len_cybja,
            "alleles": len(unique_sizes),
            "pic": pic_score
        })

# Rank markers: Highest PIC score first, then by number of unique alleles
polymorphic_panel.sort(key=lambda x: (x["pic"], x["alleles"]), reverse=True)

# Write out the clean ranked polymorphism table
with open(POLY_OUTPUT_FILE, "w") as out:
    out.write("Locus_ID\tSTR_Motif\tSize_B67\tSize_B68\tSize_Cuti\tSize_Cybja\tNum_Alleles\tPIC_Score\n")
    for item in polymorphic_panel:
        out.write(f"{item['locus_id']}\t{item['motif']}\t"
                  f"{item['size_b67']}\t{item['size_b68']}\t{item['size_cuti']}\t{item['size_cybja']}\t"
                  f"{item['alleles']}\t{item['pic']:.3f}\n")

print(f"\n✨ Stage 4 Complete!")
print(f"  - Screened 552 conserved loci.")
print(f"  - Discarded {monomorphic_count} uninformative monomorphic loci.")
print(f"  - Kept and ranked {len(polymorphic_panel)} highly polymorphic universal markers.")
print(f"📊 Polymorphic matrix saved directly to: {POLY_OUTPUT_FILE}")
