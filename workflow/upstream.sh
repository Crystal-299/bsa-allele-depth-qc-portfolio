#!/usr/bin/env bash
# Example paired-end DNA-seq path: fastp -> BWA -> samtools -> GATK gVCF -> joint SNP VCF.
# Run only on Linux with the named tools installed. No real data is bundled.
set -euo pipefail
if [[ $# -ne 3 ]]; then echo "usage: bash workflow/upstream.sh reference.fa samples.tsv output-dir" >&2; exit 2; fi
ref=$1; sheet=$2; out=$3
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 "$root/scripts/validate_upstream.py" "$ref" "$sheet"
for tool in fastp bwa samtools gatk; do command -v "$tool" >/dev/null || { echo "missing tool: $tool" >&2; exit 2; }; done
mkdir -p "$out"/{qc,bam,gvcf,variants}
samtools faidx "$ref"
dict="${ref%.*}.dict"
[[ -s "$dict" ]] || gatk CreateSequenceDictionary -R "$ref" -O "$dict"
[[ -s "${ref}.bwt" ]] || bwa index "$ref"
read -r header < "$sheet"
[[ "$header" == $'sample\tbulk\tindividuals\tploidy\tr1\tr2' ]] || { echo "invalid TSV header" >&2; exit 2; }
variants=(); high=0; low=0
while IFS=$'\t' read -r sample bulk individuals ploidy r1 r2; do
  [[ -z "$sample" ]] && continue
  [[ "$sample" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] || { echo "invalid sample: $sample" >&2; exit 2; }
  [[ "$bulk" == high || "$bulk" == low ]] || { echo "bulk must be high/low" >&2; exit 2; }
  [[ "$individuals" =~ ^[1-9][0-9]*$ && "$ploidy" =~ ^[1-9][0-9]*$ ]] || { echo "invalid pool size/ploidy" >&2; exit 2; }
  [[ -s "$r1" && -s "$r2" ]] || { echo "missing FASTQ for $sample" >&2; exit 2; }
  [[ "$bulk" == high ]] && high=$((high+1)) || low=$((low+1))
  fastp -i "$r1" -I "$r2" -o "$out/qc/${sample}_R1.fq.gz" -O "$out/qc/${sample}_R2.fq.gz" \
    --json "$out/qc/${sample}.json" --html "$out/qc/${sample}.html"
  bwa mem -t 2 -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA" "$ref" \
    "$out/qc/${sample}_R1.fq.gz" "$out/qc/${sample}_R2.fq.gz" | \
    samtools sort -@ 2 -o "$out/bam/${sample}.sorted.bam" -
  gatk MarkDuplicates -I "$out/bam/${sample}.sorted.bam" -O "$out/bam/${sample}.dedup.bam" \
    -M "$out/bam/${sample}.metrics.txt"
  samtools index "$out/bam/${sample}.dedup.bam"
  samtools flagstat -@ 2 "$out/bam/${sample}.dedup.bam" > "$out/qc/${sample}.flagstat.txt"
  samtools coverage "$out/bam/${sample}.dedup.bam" > "$out/qc/${sample}.coverage.tsv"
  gatk HaplotypeCaller -R "$ref" -I "$out/bam/${sample}.dedup.bam" \
    -O "$out/gvcf/${sample}.g.vcf.gz" -ERC GVCF --sample-ploidy "$((individuals*ploidy))"
  variants+=(--variant "$out/gvcf/${sample}.g.vcf.gz")
done < <(tail -n +2 "$sheet")
[[ $high -eq 1 && $low -eq 1 ]] || { echo "this example requires exactly one high and one low bulk" >&2; exit 2; }
gatk CombineGVCFs -R "$ref" "${variants[@]}" -O "$out/variants/combined.g.vcf.gz"
gatk GenotypeGVCFs -R "$ref" -V "$out/variants/combined.g.vcf.gz" -O "$out/variants/joint.vcf.gz"
gatk SelectVariants -R "$ref" -V "$out/variants/joint.vcf.gz" --select-type-to-include SNP \
  -O "$out/variants/snps.vcf.gz"
echo "SNP VCF: $out/variants/snps.vcf.gz (review site/genotype filters before interpretation)"
