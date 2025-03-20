#!/usr/bin/env bash
set -euo pipefail

python 0_preprocess_profiles.py
python 1_calculate_map.py
python 2_plot_map.py