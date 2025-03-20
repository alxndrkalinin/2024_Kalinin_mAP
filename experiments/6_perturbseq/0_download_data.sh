#!/bin/bash
set -euo pipefail

GSE_ID="GSE132080"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/${GSE_ID}/suppl"
OUTPUT_DIR="inputs"
mkdir -p "$OUTPUT_DIR"

# Files to download (non-"other")
declare -A FILES=(
  ["barcodes"]="${GSE_ID}_10X_barcodes.tsv.gz"
  ["data"]="${GSE_ID}_10X_matrix.mtx.gz"
  ["genes"]="${GSE_ID}_10X_genes.tsv.gz"
)
declare -A RENAMES=(
  ["barcodes"]="${GSE_ID}/barcodes.tsv"
  ["data"]="${GSE_ID}/matrix.mtx"
  ["genes"]="${GSE_ID}/genes.tsv"
)

echo "Downloading non-other files..."
for key in "${!FILES[@]}"; do
  file="${FILES[$key]}"
  output_path="${OUTPUT_DIR}/${file}"

  if [[ ! -f "${output_path}" ]]; then
    echo "Downloading ${file}..."
    wget -q -nc -P "$OUTPUT_DIR" "${BASE_URL}/${file}"
  else
    echo "File ${file} already exists, skipping."
  fi
done

# Files under "other"
OTHER_FILES=(
  "${GSE_ID}_cell_identities.csv.gz"
  "${GSE_ID}_sgRNA_barcode_sequences_and_phenotypes.csv.gz"
)

echo "Downloading 'other' files..."
for file in "${OTHER_FILES[@]}"; do
  output_path="${OUTPUT_DIR}/${file}"

  if [[ ! -f "${output_path}" ]]; then
    echo "Downloading ${file}..."
    wget -q -nc -P "$OUTPUT_DIR" "${BASE_URL}/${file}"
  else
    echo "File ${file} already exists, skipping."
  fi
done

echo "Extracting gzip files (non-other)..."
for key in "${!FILES[@]}"; do
  file="${FILES[$key]}"
  gz_path="${OUTPUT_DIR}/${file}"
  target="${OUTPUT_DIR}/$(dirname "${RENAMES[$key]}")/$(basename "${RENAMES[$key]}")"

  if [[ ! -f "${target}" ]]; then
    mkdir -p "$(dirname "$target")"
    echo "Extracting ${gz_path} to ${target}..."
    gunzip -c "$gz_path" > "$target"
  else
    echo "Extracted file ${target} already exists, skipping."
  fi
done

# Download and extract the paper supplement
PAPER_SUPP_BASE_URL="https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0387-5/MediaObjects"
PAPER_SUPP_FILENAME="41587_2019_387_MOESM3_ESM.zip"
PAPER_SUPP_DIR="paper_supplement"
PAPER_SUPP_PATH="${OUTPUT_DIR}/${PAPER_SUPP_FILENAME}"

echo "Downloading paper supplement..."
if [[ ! -f "${PAPER_SUPP_PATH}" ]]; then
  wget -q -nc -P "$OUTPUT_DIR" "${PAPER_SUPP_BASE_URL}/${PAPER_SUPP_FILENAME}"
else
  echo "File ${PAPER_SUPP_FILENAME} already exists, skipping."
fi

echo "Extracting ${PAPER_SUPP_FILENAME} to ${OUTPUT_DIR}/${PAPER_SUPP_DIR}..."
mkdir -p "${OUTPUT_DIR}/${PAPER_SUPP_DIR}"
unzip -o -q "${PAPER_SUPP_PATH}" -d "${OUTPUT_DIR}/${PAPER_SUPP_DIR}"

echo "Perturb-seq data download and extraction completed."
