#!/usr/bin/env bash
set -euo pipefail

experiments=("2_cellhealth" "3_cpg0004" "4_cpg0016orf" "5_nelisa" "6_perturbseq" "7_mitocheck")

for exp in "${experiments[@]}"; do
    echo "Processing $exp..."
    cd "experiments/$exp/"
    script_name="run_all_${exp#*_}.sh"
    bash "$script_name"
    cd ../..
done

echo "Done."
