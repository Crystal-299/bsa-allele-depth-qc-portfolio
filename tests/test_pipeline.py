import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PipelineTest(unittest.TestCase):
    def test_toy_vcf_to_plots(self):
        scratch = ROOT / "demo-output"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as folder:
            out = Path(folder)
            steps = [
                ("vcf_to_ad.py", [ROOT / "examples/toy.vcf", ROOT / "examples/project.json", out / "ad.csv"]),
                ("bsa_qc.py", [out / "ad.csv", out / "sites.csv", "--min-depth", "10", "--max-depth", "100"]),
                ("summarize.py", [out / "sites.csv", ROOT / "examples/project.json", out / "plots"]),
                ("to_qtlseqr.py", [out / "sites.csv", out / "qtlseqr.tsv"]),
            ]
            for script, args in steps:
                subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)], check=True, capture_output=True)
            with (out / "sites.csv").open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 6)
            self.assertEqual(sum(r["status"] == "pass" for r in rows), 5)
            self.assertEqual(rows[0]["delta_snp_index"], "0.500000")
            self.assertEqual(rows[0]["high_ad"], "5,15")
            for name in ("windows.csv", "delta.svg", "window_delta.svg", "depth.svg"):
                self.assertGreater((out / "plots" / name).stat().st_size, 0)
            with (out / "qtlseqr.tsv").open(encoding="utf-8") as stream:
                self.assertEqual(sum(1 for _ in stream), 6)


if __name__ == "__main__":
    unittest.main()
