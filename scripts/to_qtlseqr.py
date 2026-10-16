#!/usr/bin/env python3
"""Convert passed site AD records into QTLseqr importFromTable columns."""
import argparse
import csv
from pathlib import Path


def run(source: Path, target: Path):
    with source.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"chrom", "position", "ref", "alt", "high_ad", "low_ad", "status"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("input requires AD and status columns")
        rows = [r for r in reader if r["status"] == "pass"]
    if not rows:
        raise ValueError("no passed SNP sites")
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("CHROM", "POS", "REF", "ALT", "AD_REF.High", "AD_ALT.High", "AD_REF.Low", "AD_ALT.Low"))
        for r in rows:
            hr, ha = r["high_ad"].split(",")
            lr, la = r["low_ad"].split(",")
            writer.writerow((r["chrom"], r["position"], r["ref"], r["alt"], hr, ha, lr, la))
    print(f"sites={len(rows)} output={target}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    try:
        run(a.input, a.output)
    except (OSError, ValueError) as e:
        p.exit(2, f"error: {e}\n")


if __name__ == "__main__":
    main()
