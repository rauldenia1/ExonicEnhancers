#!/usr/bin/env bash
# ==============================================================================
# Pipeline: Identification of Exonic Enhancers (EEs) Across Species
# Description: Identifies non-redundant coding exons overlapping ATAC peaks/summits,
#              excluding promoter- and terminus-proximal regions (+/- 100 bp).
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 0. Configuration & Species Parameters
# ------------------------------------------------------------------------------
# Usage: ./identify_exonic_enhancers.sh <genome_fa> <gff_file> <atac_peaks> <chr_regex> <output_dir>
GENOME_FA="${1:-"../00_genome/MpTak1_v7.1.fa"}"
GFF_FILE="${2:-"../00_genome/MpTak_v7.1.gff"}"
ATAC_PEAKS="${3:-"../03_peaks/Marchantia_ATAC_peaks.narrowPeak"}"
CHR_REGEX="${4:-"^chr[0-9]+$"}"               # Regex for primary chromosomes
OUTPUT_DIR="${5:-"results_ee_identification"}"

GENOME_FAI="${GENOME_FA}.fai"
mkdir -p "${OUTPUT_DIR}"

echo "=== Running EE Identification Pipeline ==="
echo "Genome FASTA    : ${GENOME_FA}"
echo "GFF Annotation  : ${GFF_FILE}"
echo "ATAC Peaks      : ${ATAC_PEAKS}"
echo "Chromosome Regex: ${CHR_REGEX}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "==========================================="

# Index genome if missing
if [[ ! -f "${GENOME_FAI}" ]]; then
  echo "Generating FASTA index for ${GENOME_FA}..."
  samtools faidx "${GENOME_FA}"
fi

# ------------------------------------------------------------------------------
# 1. Reference Geometries (Filtered Chromosomes)
# ------------------------------------------------------------------------------
# Generate chromosome sizes
awk -v reg="${CHR_REGEX}" '$1 ~ reg {print $1 "\t" $2}' "${GENOME_FAI}" \
  | sort -k1,1 > "${OUTPUT_DIR}/genome.sizes"

# Extract all protein-coding genes
awk -F'\t' -v reg="${CHR_REGEX}" '
  $3 == "gene" && $1 ~ reg {
    gene_id = "."
    if (match($9, /ID=([^;]+)/, arr)) {
      gene_id = arr[1]
    }
    print $1 "\t" $4 - 1 "\t" $5 "\t" gene_id "\t.\t" $7
  }' "${GFF_FILE}" \
  | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/genes_chr.bed"

# Extract protein-coding exons and merge overlapping isoforms (strand-aware)
awk -F'\t' -v reg="${CHR_REGEX}" '$3 == "exon" && $1 ~ reg && $9 ~ /Parent=/ {print $1 "\t" $4 - 1 "\t" $5 "\t.\t.\t" $7}' "${GFF_FILE}" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -s -c 6 -o distinct -i - \
  | awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\tEE_" NR "\t.\t" $4}' > "${OUTPUT_DIR}/exons_universe.bed"

# Extract CDS features and merge overlapping isoforms (strand-aware)
awk -F'\t' -v reg="${CHR_REGEX}" '$3 == "CDS" && $1 ~ reg {print $1 "\t" $4 - 1 "\t" $5 "\t.\t.\t" $7}' "${GFF_FILE}" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -s -c 6 -o distinct -i - \
  | awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\tcds_" NR "\t.\t" $4}' > "${OUTPUT_DIR}/cds_universe.bed"

# ------------------------------------------------------------------------------
# 2. Filter Terminal Regions (TSS/TES +/- 100 bp)
# ------------------------------------------------------------------------------
awk -F'\t' -v reg="${CHR_REGEX}" '$3 == "mRNA" && $1 ~ reg {
    s = $4 - 1; 
    e = $5; 
    print $1 "\t" s "\t" s + 1 "\t.\t.\t" $7; 
    print $1 "\t" e - 1 "\t" e "\t.\t.\t" $7
  }' "${GFF_FILE}" \
  | sort -k1,1 -k2,2n \
  | bedtools slop -b 100 -g "${OUTPUT_DIR}/genome.sizes" -i - \
  | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/mrna_ends_slop100.bed"

bedtools intersect -s -v \
  -a "${OUTPUT_DIR}/exons_universe.bed" \
  -b "${OUTPUT_DIR}/mrna_ends_slop100.bed" > "${OUTPUT_DIR}/exons_filtered.bed"

# ------------------------------------------------------------------------------
# 3. ATAC-seq Peak Processing
# ------------------------------------------------------------------------------
awk -v reg="${CHR_REGEX}" '$1 ~ reg' "${ATAC_PEAKS}" \
  | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/atac_peaks.bed"

awk -F'\t' 'BEGIN {OFS="\t"} {summit = $2 + $10; print $1, summit, summit + 1, $4}' "${OUTPUT_DIR}/atac_peaks.bed" \
  > "${OUTPUT_DIR}/atac_peaks_summit.bed"

# ------------------------------------------------------------------------------
# 4. Annotation Function
# ------------------------------------------------------------------------------
annotate_ees() {
  local peak_file="$1"
  local output_prefix="$2"

  bedtools intersect -u -a "${OUTPUT_DIR}/exons_filtered.bed" -b "${peak_file}" \
    | bedtools intersect -s -a - -b "${OUTPUT_DIR}/cds_universe.bed" \
    > "${OUTPUT_DIR}/${output_prefix}.bed"

  bedtools intersect -s -wa -wb -a "${OUTPUT_DIR}/${output_prefix}.bed" -b "${OUTPUT_DIR}/genes_chr.bed" \
    | awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' \
    | bedtools groupby -g 1,2,3,4,5 -c 6 -o distinct \
    > "${OUTPUT_DIR}/${output_prefix}.annotated.bed"

  cut -f6 "${OUTPUT_DIR}/${output_prefix}.annotated.bed" \
    | tr ',' '\n' \
    | sort -u > "${OUTPUT_DIR}/target_genes_${output_prefix}.txt"
}

# ------------------------------------------------------------------------------
# 5. Execution
# ------------------------------------------------------------------------------
echo "Processing full peak overlap (Permissive)..."
annotate_ees "${OUTPUT_DIR}/atac_peaks.bed" "ee_cds_atac"

echo "Processing summit peak overlap (Stringent)..."
annotate_ees "${OUTPUT_DIR}/atac_peaks_summit.bed" "ee_cds_atac_summit"

echo "Pipeline finished successfully. Output saved in: ${OUTPUT_DIR}/"