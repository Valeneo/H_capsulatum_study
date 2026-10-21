"""Collect every sample's fastp JSON into one table: results/01_qc/qc_summary.tsv"""
import csv, json, pathlib

root = pathlib.Path(__file__).resolve().parent.parent
rows = []
for jf in sorted((root / "results/01_qc").glob("*/*.fastp.json")):
    j = json.loads(jf.read_text())
    b, a = j["summary"]["before_filtering"], j["summary"]["after_filtering"]
    rows.append({
        "run": jf.parent.name,
        "reads_raw": b["total_reads"], "reads_clean": a["total_reads"],
        "pct_reads_kept": round(100 * a["total_reads"] / b["total_reads"], 1),
        "bases_clean_mb": round(a["total_bases"] / 1e6),
        "q30_raw": round(100 * b["q30_rate"], 1), "q30_clean": round(100 * a["q30_rate"], 1),
        "gc_clean": round(100 * a["gc_content"], 1),
        "dup_rate": round(100 * j["duplication"]["rate"], 1),
        "adapter_trimmed_reads": j["adapter_cutting"].get("adapter_trimmed_reads", 0),
    })

out = root / "results/01_qc/qc_summary.tsv"
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"{len(rows)} samples -> {out}")
