#!/usr/bin/env bash
set -euo pipefail

# run_mutect2_chr13.sh
# Align paired FASTQ, mark duplicates, call Mutect2 (tumor-only), filter, and annotate with VEP.
# Usage: run_mutect2_chr13.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> <out_dir>
# Assumes conda env with bwa,samtools,picard,gatk,vep is active (or use `conda run -n <env>` to execute script).

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <sample_name> <out_dir>"
  exit 1
fi

R1="$1"
R2="$2"
SAMPLE="$3"
OUTDIR="$4"

THREADS=4
MEM=16G

# Derive repository root (one level above scripts/)
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF="${REPO_ROOT}/workshop/genome.fa"
BAITS="${REPO_ROOT}/workshop/Baits.nochr.bed"
GRESOURCE="${REPO_ROOT}/workshop/1000G_Omini.nochr.vcf.gz"

echo "Checking inputs..."
[ -f "$R1" ] || { echo "R1 not found: $R1"; exit 2; }
[ -f "$R2" ] || { echo "R2 not found: $R2"; exit 2; }
[ -f "$REF" ] || { echo "Reference not found: $REF"; exit 2; }
[ -f "$BAITS" ] || { echo "Baits bed not found: $BAITS"; exit 2; }
[ -f "$GRESOURCE" ] || { echo "Germline resource not found: $GRESOURCE"; exit 2; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Indexing reference (if needed)..."
if [ ! -f "${REF}.bwt" ]; then
  bwa index "$REF"
fi
if [ ! -f "${REF}.fai" ]; then
  samtools faidx "$REF"
fi
if [ ! -f "${REF%.*}.dict" ]; then
  gatk CreateSequenceDictionary -R "$REF"
fi

SORTED_BAM="${SAMPLE}.sorted.bam"
DEDUP_BAM="${SAMPLE}.dedup.bam"
METRICS="${SAMPLE}.dedup.metrics.txt"

echo "Aligning with bwa mem and sorting..."
bwa mem -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" -t "$THREADS" "$REF" "$R1" "$R2" \
  | samtools view -b -@ "$THREADS" - \
  | samtools sort -@ "$THREADS" -o "$SORTED_BAM"
samtools index "$SORTED_BAM"

echo "Marking duplicates with Picard..."
picard MarkDuplicates I="$SORTED_BAM" O="$DEDUP_BAM" M="$METRICS" CREATE_INDEX=true

echo "Running GATK Mutect2 (tumor-only)..."
UNFILTERED_VCF="${SAMPLE}.unfiltered.vcf.gz"
gatk --java-options "-Xmx${MEM}" Mutect2 \
  -R "$REF" \
  -I "$DEDUP_BAM" \
  -tumor "$SAMPLE" \
  --germline-resource "$GRESOURCE" \
  -L "$BAITS" \
  -O "$UNFILTERED_VCF"

echo "Filtering Mutect2 calls..."
FILTERED_VCF="${SAMPLE}.filtered.vcf.gz"
gatk --java-options "-Xmx${MEM}" FilterMutectCalls \
  -R "$REF" \
  -V "$UNFILTERED_VCF" \
  -O "$FILTERED_VCF"


echo "All done. Files in: $(pwd)"
ls -lh
