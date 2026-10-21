"""Pick phylogenetic marker genes: BUSCO genes found complete + single-copy in nearly every isolate (and the outgroup).
Writes one protein FASTA per gene to results/07_phylo/genes/. Usage: python3 07_collect_genes.py [max_genes] [min_presence]
"""
import pathlib, sys, collections

root = pathlib.Path(__file__).resolve().parent.parent
max_genes = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
min_presence = float(sys.argv[2]) if len(sys.argv) > 2 else 0.90
out = root / "results/07_phylo/genes"; out.mkdir(parents=True, exist_ok=True)

def read_faa(path):
    seq, name = [], None
    for l in open(path):
        l = l.strip()
        if l.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = l[1:], []
        elif l: seq.append(l)
    if name: yield name, "".join(seq)

# gene id -> {isolate: protein}
genes = collections.defaultdict(dict)
isolates = []
for d in sorted((root / "results/04_assembly").glob("SRR*")):
    sc = d / "busco/run_onygenales_odb10/busco_sequences/single_copy_busco_sequences"
    if not sc.exists(): continue
    isolates.append(d.name)
    for f in sc.glob("*.faa"):
        for _, s in read_faa(f): genes[f.stem][d.name] = s
og = root / "results/07_phylo/outgroup_busco/run_onygenales_odb10/busco_sequences/single_copy_busco_sequences"
for f in og.glob("*.faa"):
    for _, s in read_faa(f): genes[f.stem]["OUTGROUP_Blastomyces"] = s

n = len(isolates)
ok = [(g, d) for g, d in genes.items() if "OUTGROUP_Blastomyces" in d and sum(k != "OUTGROUP_Blastomyces" for k in d) >= min_presence * n]
# Longest, most complete genes carry the most signal per gene.
ok.sort(key=lambda gd: (-len(gd[1]), -sum(len(s) for s in gd[1].values()) / len(gd[1])))
ok = ok[:max_genes]
for g, d in ok:
    with open(out / f"{g}.faa", "w") as fh:
        for iso, s in d.items(): fh.write(f">{iso}\n{s}\n")
print(f"{n} isolates, {len(genes)} BUSCO genes seen, {len(ok)} kept (present in >= {min_presence:.0%} and in the outgroup)")
