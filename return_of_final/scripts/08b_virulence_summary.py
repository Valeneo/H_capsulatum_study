"""Combine per-isolate BLAST hits into gene presence/absence + mating type: results/08_virulence/virulence_matrix.tsv"""
import csv, pathlib

root = pathlib.Path(__file__).resolve().parent.parent
queries = [l[1:].strip() for l in open(root / "reference/virulence_queries.faa") if l.startswith(">")]
names = [q.split("|")[0] for q in queries]
MIN_ID, MIN_COV = 70.0, 60.0          # a gene is "present" if >=70% identical over >=60% of the query
MAT_ID, MAT_COV = 90.0, 90.0          # mating-type calls are stricter: a paralog hits ~40-50%, the true locus ~95%+

rows = []
for d in sorted((root / "results/08_virulence").glob("SRR*")):
    hits = {}
    f = d / "hits.tsv"
    if not f.exists(): continue
    for q, _, ident, cov, _ in csv.reader(open(f), delimiter="\t"):
        hits[q.split("|")[0]] = (float(ident), float(cov))
    row = {"run": d.name}
    for n in names:
        ident, cov = hits.get(n, (0.0, 0.0))
        row[n] = f"{ident:.0f}" if ident >= MIN_ID and cov >= MIN_COV else "absent"
    m1 = hits.get("MAT1-1-1_alpha_box", (0, 0)); m2 = hits.get("MAT1-2-1_HMG_box", (0, 0))
    has1 = m1[0] >= MAT_ID and m1[1] >= MAT_COV; has2 = m2[0] >= MAT_ID and m2[1] >= MAT_COV
    row["mating_type"] = "MAT1-1" if has1 and not has2 else "MAT1-2" if has2 and not has1 else "both?" if has1 and has2 else "unresolved"
    rows.append(row)

out = root / "results/08_virulence/virulence_matrix.tsv"
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run"] + names + ["mating_type"], delimiter="\t"); w.writeheader(); w.writerows(rows)
from collections import Counter
print(f"{len(rows)} isolates ->", out); print("mating types:", dict(Counter(r["mating_type"] for r in rows)))
