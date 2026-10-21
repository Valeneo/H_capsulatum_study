"""Verify a finished stage against explicit thresholds. Prints one line per failing sample and a summary.
Usage: python3 scripts/99_verify.py {1|2|3|4}
"""
import csv, json, pathlib, re, subprocess, sys

root = pathlib.Path(__file__).resolve().parent.parent
samples = [r["run"] for r in csv.DictReader(open(root / "config/samples.tsv"), delimiter="\t")]
stage = sys.argv[1]
bad, ok = [], 0

def fail(run, why): bad.append(run); print(f"FAIL {run}: {why}")

for run in samples:
    if stage == "1":
        j = root / f"results/01_qc/{run}/{run}.fastp.json"
        if not j.exists(): fail(run, "missing"); continue
        s = json.loads(j.read_text())["summary"]; b, a = s["before_filtering"], s["after_filtering"]
        kept = 100 * a["total_reads"] / b["total_reads"]
        if kept < 70: fail(run, f"only {kept:.0f}% of reads kept")
        elif a["q30_rate"] < 0.90: fail(run, f"Q30 after trimming {a['q30_rate']:.2f} < 0.90")
        elif a["total_bases"] < 150e6: fail(run, f"only {a['total_bases']/1e6:.0f} Mb (<5x coverage)")
        else: ok += 1
    elif stage == "2":
        d = root / f"results/02_variants/{run}"
        if not (d / f"{run}.g.vcf.gz.tbi").exists(): fail(run, "no gVCF"); continue
        m = re.search(r"\((\d+\.\d+)%", (d / f"{run}.flagstat.txt").read_text().split("primary mapped")[1])
        cov = [l.split("\t") for l in (d / f"{run}.coverage.tsv").read_text().splitlines()[1:] if "CP069113" not in l]
        depth = sum(float(c[6]) * (int(c[2]) - int(c[1]) + 1) for c in cov) / sum(int(c[2]) - int(c[1]) + 1 for c in cov)
        pct = float(m.group(1))
        if pct < 90: fail(run, f"only {pct}% of reads mapped (wrong species, or contamination?)")
        elif depth < 5: fail(run, f"nuclear depth {depth:.1f}x < 5x")
        else: ok += 1
    elif stage == "3":
        d = root / f"results/04_assembly/{run}"
        if not (d / "busco_summary.txt").exists(): fail(run, "missing"); continue
        c = float(re.search(r"C:(\d+\.\d+)%", (d / "busco_summary.txt").read_text()).group(1))
        size = int(open(d / "assembly_stats.tsv").read().splitlines()[1].split("\t")[4])
        if c < 90: fail(run, f"BUSCO complete {c}% < 90")
        elif not 25e6 < size < 45e6: fail(run, f"assembly size {size/1e6:.1f} Mb outside 25-45 Mb")
        else: ok += 1

print(f"\nstage {stage}: {ok}/{len(samples)} pass, {len(bad)} fail")
sys.exit(1 if bad else 0)
