#!/usr/bin/env python3
"""Create deterministic synthetic paired DNA reads for an upstream toolchain smoke test."""
import gzip
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "toy"
COMP = str.maketrans("ACGT", "TGCA")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = random.Random(19)
    ref = "".join(rng.choice("ACGT") for _ in range(10000))
    (OUT / "reference.fa").write_text(">toy_chr1\n" + "\n".join(ref[i:i+80] for i in range(0, len(ref), 80)) + "\n", encoding="ascii")
    for sample, alt_fraction in (("High", 0.8), ("Low", 0.2)):
        with gzip.open(OUT / f"{sample}_R1.fq.gz", "wt") as f1, gzip.open(OUT / f"{sample}_R2.fq.gz", "wt") as f2:
            for n in range(80):
                start = rng.randint(390, 450)
                fragment = list(ref[start:start+220])
                if rng.random() < alt_fraction:
                    fragment[499-start] = next(base for base in "ACGT" if base != ref[499])
                r1 = "".join(fragment[:100])
                r2 = "".join(fragment[-100:]).translate(COMP)[::-1]
                for stream, mate, seq in ((f1, 1, r1), (f2, 2, r2)):
                    stream.write(f"@{sample}_{n}/{mate}\n{seq}\n+\n{'I'*100}\n")
    print(f"synthetic reference and paired reads: {OUT}")


if __name__ == "__main__":
    main()
