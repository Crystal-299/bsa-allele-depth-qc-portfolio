#!/usr/bin/env python3
"""Validate paired bulk allele depths and compute unsmoothed SNP-index."""
import argparse
import csv
from pathlib import Path

FIELDS = ("chrom", "position", "ref", "alt", "high_ad", "low_ad")
OUTPUT = ("chrom", "position", "ref", "alt", "high_ad", "low_ad", "high_depth", "low_depth",
          "high_snp_index", "low_snp_index", "delta_snp_index", "status")


def depth(value, row_number, field):
    parts = value.split(",")
    if len(parts) != 2 or any(not part.isdecimal() for part in parts):
        raise ValueError(f"row {row_number}: {field} must be REF,ALT nonnegative integers")
    return tuple(map(int, parts))


def run(source, destination, minimum, maximum):
    if minimum < 1 or maximum < minimum:
        raise ValueError("depth range must satisfy 1 <= min-depth <= max-depth")
    with source.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or set(reader.fieldnames) != set(FIELDS):
            raise ValueError("input columns must be: " + ", ".join(FIELDS))
        seen = set()
        rows = []
        for line, record in enumerate(reader, start=2):
            if any(not record[field] for field in FIELDS):
                raise ValueError(f"row {line}: missing field")
            try:
                position = int(record["position"])
            except ValueError as exc:
                raise ValueError(f"row {line}: position must be a positive integer") from exc
            if position < 1 or str(position) != record["position"]:
                raise ValueError(f"row {line}: position must be a positive integer")
            key = (record["chrom"], position)
            if key in seen:
                raise ValueError(f"row {line}: duplicate chrom/position")
            seen.add(key)
            h_ref, h_alt = depth(record["high_ad"], line, "high_ad")
            l_ref, l_alt = depth(record["low_ad"], line, "low_ad")
            h_total, l_total = h_ref + h_alt, l_ref + l_alt
            status = "pass" if minimum <= h_total <= maximum and minimum <= l_total <= maximum else "depth_fail"
            h_index = h_alt / h_total if h_total else None
            l_index = l_alt / l_total if l_total else None
            rows.append({"chrom": record["chrom"], "position": position,
                         "ref": record["ref"], "alt": record["alt"],
                         "high_ad": record["high_ad"], "low_ad": record["low_ad"],
                         "high_depth": h_total, "low_depth": l_total,
                         "high_snp_index": f"{h_index:.6f}" if h_index is not None else "",
                         "low_snp_index": f"{l_index:.6f}" if l_index is not None else "",
                         "delta_snp_index": f"{h_index - l_index:.6f}" if h_index is not None and l_index is not None else "",
                         "status": status})
    if not rows:
        raise ValueError("input contains no sites")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} sites; {sum(row['status'] == 'pass' for row in rows)} passed depth QC; output={destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=500)
    args = parser.parse_args()
    try:
        run(args.input, args.output, args.min_depth, args.max_depth)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
