#!/bin/bash
set -euo pipefail

# Define URLs and output paths
TRAIN_URL="https://github.com/wayscience/mitocheck_data/raw/main/3.normalize_data/normalized_data/training_data__no_ic.csv.gz"
CONTROL_URL="https://zenodo.org/records/7967386/files/3.normalize_data__normalized_data.zip?download=1"
RAW_DIR="inputs"
TRAIN_OUTNAME=$(basename "$TRAIN_URL")
CONTROL_OUTNAME="3.normalized_data.zip"
TRAIN_OUTPATH="${RAW_DIR}/${TRAIN_OUTNAME}"
CONTROL_OUTPATH="${RAW_DIR}/${CONTROL_OUTNAME}"
CONTROL_UNZIP_DIR="${RAW_DIR}/normalized_data"

# Ensure the directory exists
mkdir -p "$RAW_DIR"

# Download training dataset if it doesn't exist
if [[ ! -f "$TRAIN_OUTPATH" ]]; then
    echo "Downloading training dataset: $TRAIN_OUTNAME..."
    wget -q --show-progress -O "$TRAIN_OUTPATH" "$TRAIN_URL"
else
    echo "Training dataset already exists, skipping: $TRAIN_OUTPATH"
fi

# Download control dataset if it doesn't exist
if [[ ! -f "$CONTROL_OUTPATH" ]]; then
    echo "Downloading control dataset: $CONTROL_OUTNAME..."
    wget -q --show-progress -O "$CONTROL_OUTPATH" "$CONTROL_URL"
else
    echo "Control dataset already exists, skipping: $CONTROL_OUTPATH"
fi

# Unzip control dataset only if necessary
if [[ ! -d "$CONTROL_UNZIP_DIR" ]]; then
    echo "Unzipping control dataset..."
    unzip -o -q "$CONTROL_OUTPATH" -d "$RAW_DIR"
else
    echo "Control dataset already unzipped, skipping extraction."
fi

echo "Mitocheck data download and extraction complete."
