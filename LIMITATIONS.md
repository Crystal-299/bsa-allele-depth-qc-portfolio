# Limitations / 局限

- All example variants are synthetic. No real locus or result is disclosed.
- The Python path parses VCF AD, checks depth and plots descriptive Δ windows; it does not establish a QTL or calculate G′/confidence intervals. The optional QTLseqr R script is not locally verified and requires a suitable dense real dataset and a confirmed F2 design.
- The new Linux fastp/BWA/GATK upstream script has input checks, but was not executed here because the external tools are unavailable. Synthetic FASTQs only smoke-test wiring and need not yield a meaningful variant call.
- GATK pooled ploidy is set from individuals × organism ploidy; high pooled ploidy and model assumptions require project-specific review. No generic variant hard filters or BQSR known-sites resource are supplied.
- Source BSA result files exist, but exact script-to-result provenance and data ownership need review before any original material is shared.
