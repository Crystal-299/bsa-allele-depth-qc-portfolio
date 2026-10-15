# BSA allele-depth QC and SNP-index toy

This small Python command-line tool reads per-site allele depths for two bulks, applies a declared depth range, and calculates each bulk's SNP-index and their difference. It makes the first audit step of a bulk-segregant analysis inspectable and reproducible.

## Provenance

This is a **new implementation** informed by older BSA R scripts that parse `AD`, filter read depth, and inspect allele frequencies. The source directory also contains an R history, filtered tables, QTL-style output, VCFs, and BAMs. The existence of those files supports prior analysis activity, but does not prove that each output came from a particular script version. No original sample names, variants, coordinates, unpublished loci, or sequencing data are included here.

## Requirements

- Python 3.9 or newer; standard library only.

## Run the synthetic example

From the repository root:

```sh
python scripts/bsa_qc.py examples/allele_depths.csv demo-output/sites.csv --min-depth 10 --max-depth 100
```

Expected: `demo-output/sites.csv` with 5 sites, 3 passing and 2 failing the depth filter. At `toy_chr1:100`, the high and low SNP-index values should be `0.75` and `0.25`, so `delta_snp_index` is `0.5`.

## Input contract

CSV columns: `chrom`, `position`, `ref`, `alt`, `high_ad`, `low_ad`. Each AD value is `REF,ALT` with nonnegative integer counts. Positions are positive integers; alleles must be nonempty; `(chrom, position)` pairs must be unique. A site passes when **both** bulk depths are within the inclusive limits. Failed sites are retained in the output with a reason. Malformed records fail the entire run with a line-numbered error.

## Scope

`delta_snp_index = ALT/(REF+ALT)` in the high bulk minus that in the low bulk. The script makes no QTL call, smoothing, confidence interval, multiple-testing correction, or biological interpretation. See [LIMITATIONS.md](LIMITATIONS.md).
