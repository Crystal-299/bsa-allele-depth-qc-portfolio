# BSA method and decision guide / BSA 方法与选择

## 1. Inputs / 输入

- `examples/samples.tsv`: exactly one high-phenotype and one low-phenotype bulk in this minimal example. Columns are `sample,bulk,individuals,ploidy,r1,r2`. `individuals × ploidy` is passed to GATK HaplotypeCaller as pooled sample ploidy; this follows [GATK's option definition](https://gatk.broadinstitute.org/hc/en-us/articles/360036889312-HaplotypeCaller), but large pooled ploidy has practical genotyping limits and must be reviewed.
- `examples/project.json`: VCF bulk sample names, high/low phenotype labels, population (`F2`, `BC1`, `RIL`, `F1`, or `unknown`), inheritance (`dominant`, `recessive`, `unknown`), trait-allele phase (`high_ALT`, `high_REF`, `unknown`), depth bounds and window width. These fields make assumptions **visible**, not automatically correct. Population/inheritance are not used to invent a causal locus.
- Joint VCF: `#CHROM` header with the named bulk samples and per-sample `FORMAT/AD=REF,ALT`. By default only `PASS` biallelic A/C/G/T SNPs with both ADs are used. `--allow-unfiltered` can include `FILTER=.` sites from a callset to which no site filter was applied; it does **not** certify quality. Multiallelics, indels, missing AD, and failing filter calls are skipped. For real data, inspect genotype quality, mapping quality, strand bias, duplication and missingness **before** deciding site filters.

## 2. Raw reads to VCF / 原始数据到变异

`workflow/upstream.sh` is a new, documented implementation: paired FASTQ → fastp HTML/JSON QC → BWA-MEM alignment → sorted BAM → GATK MarkDuplicates → samtools flagstat/coverage reports → per-pool HaplotypeCaller `-ERC GVCF` → CombineGVCFs → GenotypeGVCFs → SNP selection. This sequence follows the [GATK gVCF workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode). It is **not** proven to have been run in the historical BSA project, and no blanket hard-filter thresholds are baked in. A known-sites BQSR stage is omitted because an authorized known-variants resource was not provided. The shell script has not been end-to-end executed in this Windows environment.

The bundled FASTQs are deliberately tiny synthetic smoke inputs. They test wiring, not discovery power or realistic GATK genotype calls. Review `fastp` quality/adaptor reports, alignment and duplication rates, coverage/AD distributions, variant counts and low-quality flags on real data. Compare high/low pools and any parents on the **same reference build**.

## 3. Definitions / 公式

For biallelic site `i` and bulk `b`, `depth_b = AD_REF,b + AD_ALT,b` and `SNPindex_b = AD_ALT,b / depth_b`. The sign convention is `ΔSNP-index = SNPindex_high − SNPindex_low`. A fixed genomic window reports the unweighted mean Δ of passed sites; this is a descriptive plot, **not** a confidence interval or QTL call. Both bulks must pass configured inclusive depth bounds. These definitions match the [QTLseqr vignette](https://github.com/bmansfeld/QTLseqr/blob/master/vignettes/QTLseqr.Rmd), but this Python summary is intentionally simpler than QTLseqr's smoothing and simulations.

## 4. Which statistic? / 指标选择

| Question / 场景 | Use / 选择 | Caution / 注意 |
| --- | --- | --- |
| Are AD values and depth credible? / 先看位点质量 | Depth/AD QC and per-site SNP-index | Histograms and hard filters depend on library and genome; no universal 10–500 threshold. |
| Which direction separates two extreme bulks? / 两池差异方向 | Signed ΔSNP-index, plus per-bulk values | REF/ALT are relative to the reference, not automatically the causal parent. |
| Is there a broad genomic signal? / 连续区段趋势 | Window means or QTLseqr smoothed Δ | Sparse SNPs and varying density can make windows misleading. |
| Is a dense-SNP signal statistically unusual? / 显著性 | QTLseqr G′/q-values, and for a supported F2 design its simulated Δ intervals | Requires appropriate population, bulk sizes, window width and enough SNPs; tiny fixture is unsuitable. |

QTLseqr's `importFromTable`, `runGprimeAnalysis`, `runQTLseqAnalysis` and plotting API are documented in its [official vignette](https://github.com/bmansfeld/QTLseqr/blob/master/vignettes/QTLseqr.Rmd). The optional R script is limited to a confirmed F2 population; its statistical branch was not tested here. Do not apply F2 simulation intervals to an F1/BC1/RIL experiment simply because a dropdown allows those labels.

For a real, sufficiently dense F2 data set, install QTLseqr following its vignette, export the passed sites with `scripts/to_qtlseqr.py`, then run `Rscript scripts/qtlseqr_optional.R <dense_sites.tsv> F2 <high_bulk_individuals> <low_bulk_individuals> <window_bp> <output_dir>`. The tiny bundled fixture intentionally fails the minimum-SNP guard and must not be used to make significance claims. The R branch was not executed in this environment.

## 5. Genetic interpretation / 遗传背景与显隐性

Record cross direction, parents, ploidy, segregating generation, bulk selection rule, affected/unaffected labels and bulk individual counts **before** interpreting sign or expected SNP-index. In a clean recessive F2 model, affected plants may be enriched for homozygous causal alleles, whereas a dominant affected bulk can mix heterozygotes and homozygotes. Those are design hypotheses, not universal numerical filters; linkage, incomplete penetrance, sampling and reference allele orientation can reverse or dilute a signal. Without parent genotype/phase, `trait_allele_phase` should remain `unknown`. Neither Δ nor G′ alone establishes causality—follow up candidate variants with segregation and independent validation.

## 6. Provenance / 来源

The historical local BSA directory contains `AD.R`, `BSA_downstream.R`, `QTLseqr.R`, VCF/BAM and result files. `AD.R` explores AD-based frequency and filtering; `QTLseqr.R` calls QTLseqr G′/Δ analyses. The exact code-to-result chain and upstream raw-read commands were not established. This public repository's upstream shell, VCF parser, plots, fixtures and validations are **new implementations** inspired by that evidence, not a republication of the original data or proof of a historical end-to-end run.
