#!/bin/bash
set -euo pipefail

bash 0_download_data.sh
python 1_preprocess_mitocheck.py
python 2_calculate_map.py
python 3_plot_map.py
