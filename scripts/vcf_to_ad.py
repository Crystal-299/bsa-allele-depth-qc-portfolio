#!/usr/bin/env python3
"""Extract biallelic SNP AD for two named bulks from a VCF/VCF.gz."""
import argparse
import csv
import gzip
import json
from pathlib import Path

FIELDS = ("chrom", "position", "ref", "alt", "high_ad", "low_ad")


def run(vcf: Path, config: Path, output: Path, allow_unfiltered=False):
    cfg = json.loads(config.read_text(encoding="utf-8"))
    high, low = cfg["high_sample"], cfg["low_sample"]
    if not high or not low or high == low:
        raise ValueError("high_sample and low_sample must name different VCF samples")
    opener = gzip.open if vcf.suffix == ".gz" else open
    rows, skipped, unfiltered, header = [], 0, 0, None
    with opener(vcf, "rt", encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n\r").split("\t")
                if high not in header[9:] or low not in header[9:]:
                    raise ValueError("configured bulk names not found in VCF sample header")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError("VCF #CHROM header missing")
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(header):
                raise ValueError(f"line {line_no}: VCF column count mismatch")
            chrom, pos, _, ref, alt, _, filt, _, fmt = parts[:9]
            if (filt != "PASS" and not (allow_unfiltered and filt == ".")) or len(ref) != 1 or len(alt) != 1 or alt not in "ACGT" or ref not in "ACGT":
                skipped += 1
                continue
            if filt == ".":
                unfiltered += 1
            keys = fmt.split(":")
            if "AD" not in keys:
                skipped += 1
                continue
            ad_index = keys.index("AD")
            def ad(sample):
                values = parts[header.index(sample)].split(":")
                value = values[ad_index] if ad_index < len(values) else "."
                depths = value.split(",")
                if len(depths) != 2 or any(not item.isdecimal() for item in depths):
                    return None
                return value
            high_ad, low_ad = ad(high), ad(low)
            if high_ad is None or low_ad is None:
                skipped += 1
                continue
            position = int(pos)
            if position < 1:
                raise ValueError(f"line {line_no}: invalid position")
            rows.append(dict(zip(FIELDS, (chrom, position, ref, alt, high_ad, low_ad))))
    if header is None or not rows:
        raise ValueError("no usable PASS biallelic SNP records with two AD values")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"retained={len(rows)} skipped={skipped} unfiltered_included={unfiltered} output={output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--allow-unfiltered", action="store_true", help="include FILTER=. sites from an unfiltered callset; review QC before interpretation")
    args = parser.parse_args()
    try:
        run(args.vcf, args.config, args.output, args.allow_unfiltered)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
