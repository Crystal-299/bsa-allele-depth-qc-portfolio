#!/usr/bin/env python3
"""Validate reference and pooled FASTQ sample sheet before expensive tools run."""
import csv
import re
import sys
from pathlib import Path


def validate(reference: Path, sheet: Path):
    if not reference.is_file() or reference.stat().st_size == 0:
        raise ValueError("reference FASTA missing/empty")
    with sheet.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = ("sample", "bulk", "individuals", "ploidy", "r1", "r2")
        if tuple(reader.fieldnames or ()) != fields:
            raise ValueError("sample sheet columns must be: " + ", ".join(fields))
        rows = list(reader)
    if len(rows) != 2 or {r["bulk"] for r in rows} != {"high", "low"}:
        raise ValueError("this minimal workflow requires exactly one high and one low pool")
    if len({r["sample"] for r in rows}) != 2:
        raise ValueError("sample IDs must be unique")
    for r in rows:
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", r["sample"] or ""):
            raise ValueError("sample IDs must be simple ASCII tokens")
        for key in ("individuals", "ploidy"):
            if not (r[key] or "").isdecimal() or int(r[key]) < 1:
                raise ValueError(f"{key} must be a positive integer")
        for key in ("r1", "r2"):
            path = Path(r[key])
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"missing/empty FASTQ: {path}")
    print("upstream inputs valid: 2 pools and 4 FASTQ files")


if __name__ == "__main__":
    try:
        validate(Path(sys.argv[1]), Path(sys.argv[2]))
    except (IndexError, OSError, ValueError) as exc:
        sys.exit(f"error: {exc}")
