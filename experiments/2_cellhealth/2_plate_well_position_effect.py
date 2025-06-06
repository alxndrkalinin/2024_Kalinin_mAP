#!/usr/bin/env python3
# coding: utf-8

# Assess plate and well position effects (Fig 3A)
#
# Use mAP to assess plate and well position effects by cross-checking phenotypic activity retrieval (replicates against neg cons) rates for 3 scenarios:
#  * replicates located on different plates and in different well positions
#  * replicates located on different plates, but in the same well position
#  * replicates located on the same plate and in different well positions

from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd


from cell_health_utils import (
    subset_6_replicates,
    get_cell_line_colors,
    stouffer_method,
)
from map_utils.map import calculate_map
from map_utils.plot import plot_map_scatter_kde, set_plotting_style, save_plot


def compute_map_results(data_files, cell_lines, mode_configs, map_config):
    """Compute mAP for plate/well position effects."""
    results = []

    for method, filepath in data_files.items():
        df = pd.read_csv(filepath, sep="\t")
        for cell_line in cell_lines:
            df_6wells = subset_6_replicates(df, cell_line)
            plate_pairs = list(combinations(df_6wells["Metadata_Plate"].unique(), 2))

            for plate_pair in plate_pairs:
                df_subset = df_6wells.query(
                    "Metadata_Plate in @plate_pair"
                ).reset_index(drop=True)

                for mode, config in mode_configs.items():
                    pair_config = {
                        "pos_sameby": ["Metadata_pert_name", "Metadata_control_index"],
                        "pos_diffby": [],
                        "neg_sameby": [],
                        "neg_diffby": ["Metadata_pert_name", "Metadata_control_index"],
                    }
                    pair_config["pos_sameby"].extend(config["pos_sameby"])
                    pair_config["pos_diffby"].extend(config["pos_diffby"])

                    map_results = calculate_map(df_subset, pair_config, map_config)
                    map_results["mode"] = mode
                    map_results["Preprocessing"] = method
                    map_results["Cell type"] = cell_line
                    results.append(map_results)

    return pd.concat(results).reset_index(drop=True)


def aggregate_results(results_df, results_file):
    """Aggregate results and save to CSV."""
    results_agg = (
        results_df.groupby(["Metadata_pert_name", "mode", "Preprocessing", "Cell type"])
        .agg({"mAP": "mean", "p_value": stouffer_method})
        .reset_index()
    )
    results_agg["p < 0.05"] = results_agg["p_value"] < 0.05
    results_agg["markers"] = np.where(
        results_agg["p_value"] < 0.05, "p<0.05", "p>=0.05"
    )
    results_agg["-log10(mAP p-value)"] = -np.log10(
        results_agg["p_value"].clip(lower=np.finfo(float).eps)
    )
    results_file.parent.mkdir(parents=True, exist_ok=True)
    results_agg.to_csv(results_file, index=False)
    print(f"Saved aggregated results to: {results_file}")
    return results_agg


def plot_plate_well_map_scatter(results_agg, figure_dir, style_kwargs=None):
    """Plot mAP effects for plate/well position effects."""
    set_plotting_style(**style_kwargs if style_kwargs else {})
    cell_line_colors = get_cell_line_colors()
    fig = plot_map_scatter_kde(
        results_agg,
        "mode",
        "",
        hue_col="Cell type",
        palette=cell_line_colors,
        y_label="Preprocessing",
        row="Preprocessing",
        legend_loc="center left",
        pr_x=0.35,
        pr_y=0.9,
        figure="Fig3A",
        size_rescale=0.28,
        point_size=5,
        save_path=None,  # We handle saving manually
    )
    save_plot(fig, "Fig3A_plate_well_position_effects", figure_dir)


def main():
    output_dir = Path("outputs")
    figure_dir = output_dir / "figures"
    results_file = output_dir / "plate_well_position_map.csv"

    data_files = {
        "Standardize": output_dir
        / "cell_health_profiles_merged_standardized_featureselected.tsv.gz",
        "MAD robustize": output_dir
        / "cell_health_profiles_merged_wholeplate_normalized_featureselected.tsv.gz",
    }

    cell_lines = ["A549", "ES2", "HCC44"]

    mode_configs = {
        "different plate, different well": {
            "pos_sameby": [],
            "pos_diffby": ["Metadata_Plate", "Metadata_Well"],
        },
        "different plate, same well": {
            "pos_sameby": ["Metadata_Well"],
            "pos_diffby": ["Metadata_Plate"],
        },
        "same plate, different well": {
            "pos_sameby": ["Metadata_Plate"],
            "pos_diffby": ["Metadata_Well"],
        },
    }

    map_config = {
        "null_size": 1_000_000,
        "groupby_columns": ["Metadata_pert_name"],
    }

    if results_file.exists():
        print(f"Loading existing results from: {results_file}")
        results_agg = pd.read_csv(results_file)
    else:
        results_df = compute_map_results(data_files, cell_lines, mode_configs, map_config)
        results_agg = aggregate_results(results_df, results_file)

    plot_plate_well_map_scatter(results_agg, figure_dir, style_kwargs={"font_size": 5, "linewidth": 0.35})


if __name__ == "__main__":
    main()
