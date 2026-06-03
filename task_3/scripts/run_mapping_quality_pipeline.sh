#!/usr/bin/env bash

set -e

SAMPLE="SRR001666"

REF="data/reference/ecoli_reference.fna"
READ1="data/fastq/${SAMPLE}_1.fastq.gz"
READ2="data/fastq/${SAMPLE}_2.fastq.gz"

FASTQC_DIR="results/fastqc"
ALIGN_DIR="results/alignment"
FLAGSTAT_DIR="results/flagstat"
VARIANT_DIR="results/variants"
LOG_DIR="logs"

SAM_FILE="${ALIGN_DIR}/${SAMPLE}.sam"
BAM_FILE="${ALIGN_DIR}/${SAMPLE}.bam"
SORTED_BAM_FILE="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
FLAGSTAT_FILE="${FLAGSTAT_DIR}/${SAMPLE}_flagstat.txt"
MAPPING_RESULT_FILE="${FLAGSTAT_DIR}/${SAMPLE}_mapping_quality.txt"
VCF_FILE="${VARIANT_DIR}/${SAMPLE}.vcf"

mkdir -p "$FASTQC_DIR" "$ALIGN_DIR" "$FLAGSTAT_DIR" "$VARIANT_DIR" "$LOG_DIR"

echo "Step 1: Running FastQC..."
fastqc "$READ1" "$READ2" -o "$FASTQC_DIR" > "${LOG_DIR}/fastqc.log" 2>&1

echo "Step 2: Indexing reference genome..."
if [ ! -f "${REF}.bwt" ]; then
  bwa index "$REF" > "${LOG_DIR}/bwa_index.log" 2>&1
else
  echo "Reference index already exists, skipping bwa index."
fi

echo "Step 3: Running BWA MEM..."
bwa mem "$REF" "$READ1" "$READ2" > "$SAM_FILE" 2> "${LOG_DIR}/bwa_mem.log"

echo "Step 4: Converting SAM to BAM..."
samtools view -S -b "$SAM_FILE" > "$BAM_FILE"

echo "Step 5: Running samtools flagstat..."
samtools flagstat "$BAM_FILE" > "$FLAGSTAT_FILE"

echo "Step 6: Parsing mapped reads percentage..."
MAPPED_PERCENT=$(grep " mapped (" "$FLAGSTAT_FILE" | head -1 | sed -E 's/.*\(([0-9.]+)% :.*/\1/')

echo "Mapped reads: ${MAPPED_PERCENT}%" | tee "$MAPPING_RESULT_FILE"

QUALITY_STATUS=$(awk -v p="$MAPPED_PERCENT" 'BEGIN {
  if (p > 90) {
    print "OK"
  } else {
    print "not OK"
  }
}')

echo "$QUALITY_STATUS" | tee -a "$MAPPING_RESULT_FILE"

if [ "$QUALITY_STATUS" = "not OK" ]; then
  echo "Mapping quality is not OK. Pipeline stopped."
  exit 0
fi

echo "Mapping quality is OK. Continuing pipeline..."

echo "Step 7: Sorting BAM..."
samtools sort "$BAM_FILE" -o "$SORTED_BAM_FILE"

echo "Step 8: Indexing sorted BAM..."
samtools index "$SORTED_BAM_FILE"

echo "Step 9: Running FreeBayes if available..."
if command -v freebayes >/dev/null 2>&1; then
  freebayes -f "$REF" "$SORTED_BAM_FILE" > "$VCF_FILE"
  echo "Variant calling finished."
  echo "VCF file: $VCF_FILE"
else
  echo "FreeBayes is not installed. Skipping variant calling."
  echo "FreeBayes is not installed. Variant calling was skipped." > "${VARIANT_DIR}/freebayes_skipped.txt"
fi

echo "Finished"