import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PipelineTest(unittest.TestCase):
    def test_unfiltered_sites_require_opt_in(self):
        scratch = ROOT / "demo-output"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as folder:
            out = Path(folder)
            text = (ROOT / "examples/toy.vcf").read_text(encoding="utf-8")
            (out / "unfiltered.vcf").write_text(text.replace("60\tPASS\t", "60\t.\t", 1), encoding="utf-8")
            base = [sys.executable, str(ROOT / "scripts/vcf_to_ad.py"), str(out / "unfiltered.vcf"),
                    str(ROOT / "examples/project.json")]
            subprocess.run([*base, str(out / "strict.csv")], check=True, capture_output=True)
            subprocess.run([*base, str(out / "inclusive.csv"), "--allow-unfiltered"], check=True, capture_output=True)
            with (out / "strict.csv").open(newline="", encoding="utf-8") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 5)
            with (out / "inclusive.csv").open(newline="", encoding="utf-8") as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 6)

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
