#!/usr/bin/env bash
set -euo pipefail

bash 0_download_data.sh
conda run -n perturbseq-processing Rscript 1_process_perturbseq.r
python 2_finalize_perturbseq.py
python 3_calculate_map_perturbseq.py
python 4_plot_map_perturbseq.py
