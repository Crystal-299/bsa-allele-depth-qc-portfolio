#!/usr/bin/env python3
"""Aggregate passed SNP-index sites in fixed windows; render simple SVG QC plots."""
import argparse
import csv
import html
import json
from collections import defaultdict
from pathlib import Path


def save_svg(path, title, points, y_min, y_max, color):
    width, height = 800, 300
    dots = []
    for i, (label, value) in enumerate(points):
        x = 60 + i * 700 / max(1, len(points) - 1)
        y = 255 - (value - y_min) * 210 / (y_max - y_min)
        dots.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{color}"><title>{html.escape(str(label))}: {value:.3f}</title></circle>')
    body = "\n".join(dots)
    ticks = "".join(f'<text x="4" y="{255-(v-y_min)*210/(y_max-y_min)+4:.1f}" font-size="11">{v:.1f}</text>'
                    for v in (y_min, (y_min + y_max)/2, y_max))
    path.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
                    f'<rect width="800" height="300" fill="white"/><text x="55" y="27" font-size="18">{html.escape(title)}</text>'
                    '<path d="M55 45V255H770" fill="none" stroke="#333"/>' + ticks + body +
                    '<text x="300" y="285" font-size="12">Sites/windows ordered by chromosome and position</text></svg>\n', encoding="utf-8")


def run(source: Path, config: Path, outdir: Path):
    cfg = json.loads(config.read_text(encoding="utf-8"))
    window = int(cfg["window_bp"])
    if window < 1:
        raise ValueError("window_bp must be positive")
    if cfg["population"] not in ("F2", "BC1", "RIL", "F1", "unknown"):
        raise ValueError("population must be F2, BC1, RIL, F1 or unknown")
    if cfg["inheritance"] not in ("dominant", "recessive", "unknown"):
        raise ValueError("inheritance must be dominant, recessive or unknown")
    if cfg["trait_allele_phase"] not in ("high_ALT", "high_REF", "unknown"):
        raise ValueError("trait_allele_phase must be high_ALT, high_REF or unknown")
    with source.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        need = {"chrom", "position", "high_depth", "low_depth", "high_snp_index", "low_snp_index", "delta_snp_index", "status"}
        if not need.issubset(reader.fieldnames or []):
            raise ValueError("site CSV lacks required columns")
        rows = [r for r in reader if r["status"] == "pass"]
    if not rows:
        raise ValueError("no passed sites to plot")
    grouped = defaultdict(list)
    for row in rows:
        pos = int(row["position"])
        start = ((pos - 1) // window) * window + 1
        grouped[(row["chrom"], start)].append(float(row["delta_snp_index"]))
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "windows.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("chrom", "start", "end", "n_sites", "mean_delta_snp_index"))
        for (chrom, start), values in sorted(grouped.items()):
            writer.writerow((chrom, start, start + window - 1, len(values), f"{sum(values)/len(values):.6f}"))
    ordered = sorted(rows, key=lambda r: (r["chrom"], int(r["position"])))
    save_svg(outdir / "delta.svg", "Delta SNP-index (passed sites)",
             [(f'{r["chrom"]}:{r["position"]}', float(r["delta_snp_index"])) for r in ordered], -1, 1, "#18766b")
    save_svg(outdir / "window_delta.svg", "Mean Delta SNP-index by fixed window",
             [(f"{chrom}:{start}", sum(values)/len(values)) for (chrom, start), values in sorted(grouped.items())],
             -1, 1, "#714e9e")
    save_svg(outdir / "depth.svg", "Depth QC (high + low)",
             [(f'{r["chrom"]}:{r["position"]}', int(r["high_depth"]) + int(r["low_depth"])) for r in ordered],
             0, max(1, max(int(r["high_depth"]) + int(r["low_depth"]) for r in ordered)), "#a25b16")
    print(f"passed_sites={len(rows)} windows={len(grouped)} output={outdir}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sites", type=Path)
    p.add_argument("config", type=Path)
    p.add_argument("output_dir", type=Path)
    a = p.parse_args()
    try:
        run(a.sites, a.config, a.output_dir)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as e:
        p.exit(2, f"error: {e}\n")


if __name__ == "__main__":
    main()
