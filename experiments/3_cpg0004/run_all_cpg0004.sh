#!/usr/bin/env bash
set -euo pipefail

snakemake -c1
python 1_plot_map.py
