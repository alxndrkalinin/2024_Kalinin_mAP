#!/usr/bin/env bash
set -euo pipefail

commit="cd91bd0daacef2b5ea25dcceb62482bb664d9de1"
output_dir="inputs"
mkdir -p "${output_dir}"

plates=(
    SQ00014610 SQ00014611 SQ00014612
    SQ00014613 SQ00014614 SQ00014615
    SQ00014616 SQ00014617 SQ00014618
)

for plate in "${plates[@]}"; do
    url="https://github.com/broadinstitute/cell-health/raw/${commit}/1.generate-profiles/data/profiles/${plate}/${plate}_augmented.csv.gz"
    output_file="${output_dir}/${plate}_augmented.csv.gz"

    if [[ ! -f "${output_file}" ]]; then
        echo "Downloading ${plate}..."
        wget -q -O "${output_file}" "${url}"
    else
        echo "File ${output_file} already exists, skipping."
    fi
done

# CERES data
figshare_base_url="https://ndownloader.figshare.com/files"
declare -A files=(
    ["ceres.csv"]="24613292"
    ["depmap_sample_info.csv"]="24613394"
)

for file in "${!files[@]}"; do
    url="${figshare_base_url}/${files[$file]}"
    output_file="${output_dir}/${file}"

    if [[ ! -f "${output_file}" ]]; then
        echo "Downloading ${file}..."
        wget -q -O "${output_file}" "${url}"
    else
        echo "File ${output_file} already exists, skipping."
    fi
done

echo "Cell Health data download complete."