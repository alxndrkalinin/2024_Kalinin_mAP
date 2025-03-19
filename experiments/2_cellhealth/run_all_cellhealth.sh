#!/usr/bin/env bash
set -euo pipefail

bash 0_download_data.sh
python 1_preprocess_data.py
python 2_plate_well_position_effect.py
python 3_channel_compartment_analysis.py
python 4_phenotypic_activity.py
python 5_phenotypic_consistency.py
