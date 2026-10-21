#!/usr/bin/env bash
# Objective 4 (phylogenomics): BUSCO single-copy genes -> align -> trim -> concatenate -> partitioned ML tree.
# Usage: bash scripts/07_phylogeny.sh [threads] [max_genes]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$HOME/miniconda3/envs/hcap/bin:$PATH"
THREADS="${1:-20}"; MAXG="${2:-1000}"
W="$ROOT/results/07_phylo"; mkdir -p "$W"; cd "$W"

# 0. The outgroup (Blastomyces dermatitidis, the closest well-annotated relative) goes through the same BUSCO gene set,
#    giving every gene tree a root: it tells us which split is the OLDEST, so the tree has a direction.
if [[ ! -d outgroup_busco ]]; then
  busco -i "$ROOT/reference/proteomes/Blastomyces_dermatitidis_outgroup.faa" -m proteins -l onygenales_odb10 -c "$THREADS" \
    --download_path "$ROOT/busco_downloads" --offline -o outgroup_busco --out_path "$W" -f > outgroup_busco.log 2>&1
fi

# 1. Choose genes that are single-copy in >=90% of isolates and in the outgroup.
python3 "$ROOT/scripts/07_collect_genes.py" "$MAXG" 0.90

# 2. Align each gene separately (MAFFT --auto picks the method for the data size). One gene = one independent evolutionary history.
mkdir -p aln trim
ls genes/*.faa | xargs -P "$THREADS" -I{} sh -c 'b=$(basename {} .faa); [ -s aln/$b.aln ] || mafft --auto --thread 1 --quiet {} > aln/$b.aln'

# 3. Trim alignment columns that are mostly gaps or uninformative (ClipKIT smart-gap). Poorly aligned columns add noise, not signal.
ls aln/*.aln | xargs -P "$THREADS" -I{} sh -c 'b=$(basename {} .aln); [ -s trim/$b.clip ] || clipkit {} -m smart-gap -o trim/$b.clip -q'

# 4. Concatenate into one supermatrix; write a partition file so each gene keeps its OWN substitution model.
python3 - <<'PY'
import glob, os, collections
names, parts, pos = None, [], 1
seqs = collections.defaultdict(list)
files = sorted(glob.glob("trim/*.clip"))
alns = []
for f in files:
    d, n, s = {}, None, []
    for l in open(f):
        l = l.strip()
        if l.startswith(">"):
            if n: d[n] = "".join(s)
            n, s = l[1:], []
        elif l: s.append(l)
    if n: d[n] = "".join(s)
    L = len(next(iter(d.values())))
    if L < 50: continue
    alns.append((os.path.basename(f)[:-5], d, L))
taxa = sorted({t for _, d, _ in alns for t in d})
for gene, d, L in alns:
    for t in taxa: seqs[t].append(d.get(t, "-" * L))     # a gene missing in a taxon is filled with gaps
    parts.append(f"  charset {gene} = {pos}-{pos+L-1};"); pos += L
with open("supermatrix.faa", "w") as fh:
    for t in taxa: fh.write(f">{t}\n{''.join(seqs[t])}\n")
with open("partitions.nex", "w") as fh:
    fh.write("#nexus\nbegin sets;\n" + "\n".join(parts) + "\nend;\n")
print(f"supermatrix: {len(taxa)} taxa x {pos-1} amino-acid sites, {len(alns)} genes")
PY

# 5. Maximum-likelihood tree.
#    -p partitions   each gene is its own partition      -m MFP+MERGE  test substitution models and merge similar genes
#    -rcluster 10    only test the 10% most promising merges (keeps it tractable)
#    -B 1000 -bnni   1000 ultrafast bootstraps, with NNI refinement to avoid inflated support
#    -o OUTGROUP     root the tree on Blastomyces
iqtree2 -s supermatrix.faa -p partitions.nex -m MFP+MERGE -rcluster 10 -B 1000 -bnni \
  -T "$THREADS" -o OUTGROUP_Blastomyces --prefix hcap_tree --redo > iqtree.log 2>&1
tail -5 hcap_tree.log || true
