# BSA workflow portfolio | BSA 分析流程作品集

## 中文

本仓库把双端 DNA-seq 的 **fastp 质控 → BWA 比对 → GATK gVCF/联合分型 → SNP/AD 提取 → 深度过滤 → SNP-index/ΔSNP-index → 窗口图**串成可检查的流程。还提供 QTLseqr 的可选 R 入口。原 BSA 目录有 AD、QTLseqr 脚本和结果文件，但无法证明它们与某版完整上游流程逐步对应；本仓库的上游、VCF 解析、图表和玩具数据均是**本次新整理实现**，不包含真实测序数据或未发表位点。[方法和指标选择](docs/METHODS.md)。

**无需生信软件的可运行小样例（Python 3.9+）：**

```bash
python scripts/vcf_to_ad.py examples/toy.vcf examples/project.json demo-output/ad.csv
python scripts/bsa_qc.py demo-output/ad.csv demo-output/sites.csv --min-depth 10 --max-depth 100
python scripts/summarize.py demo-output/sites.csv examples/project.json demo-output/plots
python scripts/to_qtlseqr.py demo-output/sites.csv demo-output/qtlseqr.tsv
python -m unittest discover -s tests -v
```

预期：7 条虚构 VCF 记录中保留 6 条 PASS 双等位 SNP；5 条通过深度过滤；`toy_chr1:100` 的两池指数为 0.75/0.25，Δ=0.5；输出 `windows.csv`、`delta.svg`、`window_delta.svg`、`depth.svg` 和 QTLseqr 输入表。参数入口在 `examples/project.json`：两个池的样本/表型、群体背景、显隐性、等位方向、深度范围和窗口宽度。**群体/显隐性字段只记录解释前提，不自动判定因果。**

**从 FASTQ 开始（Linux；需 fastp、BWA、samtools、GATK4）：**

```bash
python scripts/make_toy_fastq.py
python scripts/validate_upstream.py tests/toy/reference.fa examples/samples.tsv
bash workflow/upstream.sh tests/toy/reference.fa examples/samples.tsv tests/toy/upstream
python scripts/vcf_to_ad.py tests/toy/upstream/variants/snps.vcf.gz examples/project.json demo-output/upstream_ad.csv
```

样本表列为 `sample, bulk, individuals, ploidy, r1, r2`；真实分析须替换成**获授权且版本匹配**的参考基因组和 FASTQ，先审查 fastp、比对、重复率及变异过滤。此 FASTQ→GATK 段在当前 Windows 环境**未实跑**，玩具数据也不保证能产生有意义的候选区段。下游小样例与输入检查已运行通过。QTLseqr 的 G′/F2 区间统计仅适用于合适的真实设计和足够 SNP，不在小样例中声称显著性。[局限](LIMITATIONS.md)。

## English

This repository connects read QC (fastp), BWA alignment, GATK gVCF/joint calling, VCF allele depths, depth gates, SNP-index/ΔSNP-index, window summaries and SVG plots. The historical BSA folder evidences AD/QTLseqr analysis, **not** a verified end-to-end upstream command chain. The public pipeline and fixtures are new reconstructions; no original reads or loci are shared.

Run the five Python commands above to reproduce the synthetic downstream example (6 retained SNPs, 5 depth-passing sites, Δ=0.5 at `toy_chr1:100`). The four Linux commands show the optional raw-read path with external tools; that stage was not executed here. Configure sample IDs, pool sizes and ploidy in `examples/samples.tsv`, and population, inheritance, allele phase and depth/window choices in `examples/project.json`. See [methods](docs/METHODS.md) for formulas, metric choice and genetic interpretation; see [limitations](LIMITATIONS.md) for unverified stages.
