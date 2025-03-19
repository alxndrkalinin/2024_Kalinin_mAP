#!/usr/bin/env bash
set -euo pipefail

bash 0_download_data.sh
python 1_phenotypic_activity_orf.py 
python 2_phenotypic_consistency_corum_complex.py
python 3_plot_map_results.py