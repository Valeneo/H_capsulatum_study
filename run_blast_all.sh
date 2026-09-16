#!/bin/bash
set -euo pipefail

# Directory paths
GENOMES_DIR="genomes"
FLANK_DIR="flanking_sequences"
BLAST_OUT_DIR="blast/outputs"
DB_DIR="blast/databases"
THREADS=14

mkdir -p "${BLAST_OUT_DIR}"
mkdir -p "${DB_DIR}"

echo "=========================================================================="
echo "⚡ Histoplasma capsulatum Cross-Genome BLASTn Search (${THREADS} Threads)"
echo "=========================================================================="

QUERY_FASTA="${FLANK_DIR}/Hcaps_G217B_Reference_consensus_flanks.fasta"

if [ ! -f "${QUERY_FASTA}" ]; then
    echo "❌ Error: Query flank FASTA not found at ${QUERY_FASTA}."
    exit 1
fi

echo "🔍 Query Dataset: $(basename "${QUERY_FASTA}")"

# 1. Build BLAST databases for all target genomes
for genome in "${GENOMES_DIR}"/*.fasta; do
    base=$(basename "${genome}" .fasta)
    db_prefix="${DB_DIR}/${base}"
    
    echo "📦 [1/2] Creating BLAST database for: ${base}..."
    makeblastdb -in "${genome}" -dbtype nucl -out "${db_prefix}" -title "${base}" > /dev/null
    
    # 2. Run BLASTn of query flanks against each target database
    if [ "${base}" != "Hcaps_G217B_Reference" ]; then
        tab_out="${BLAST_OUT_DIR}/G217B_vs_${base}.tab"
        echo "🚀 [2/2] Running BLASTn: G217B vs ${base} (14 threads)..."
        blastn \
            -query "${QUERY_FASTA}" \
            -db "${db_prefix}" \
            -out "${tab_out}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -evalue 1e-5 \
            -max_target_seqs 5 \
            -num_threads "${THREADS}"
        echo "   ✅ Output saved: ${tab_out} ($(wc -l < "${tab_out}") hits)"
    fi
done

echo "=========================================================================="
echo "🎉 BLASTn SEARCH FINISHED! You are ready for Step 3."
echo "=========================================================================="

