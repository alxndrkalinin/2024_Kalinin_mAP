#!/usr/bin/env python3
# coding: utf-8

# 1. Plot mAP vs alternative metrics calculated from simulated data (Figures 2, S2-5)
#
#
# We generated perturbation and control profiles such that each perturbation had 2 to 4 replicates,
# in experiments with 12, 24, or 36 control replicates (in all combinations).
# The simulated profiles varied in feature size, ranging from 100 to 5,000, representing a typical range
# for various kinds of profiles, particularly after feature reduction or selection.
# Specificed number of features between 1 and 100% were perturbed (i.e. sampled from a different distribution).

from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from map_utils.plot import set_plotting_style, save_plot
from simulation_utils import (
    get_figure2_palette,
    get_figureS5_palette,
    plot_simulation_results,
)


def plot_and_save(
    csv_path, metrics, palette, figure_name, output_dir, plot_params=None
):
    """
    Load simulated data from CSV, plot simulation results, and save the figure.
    """
    data = pd.read_csv(csv_path)
    plot_params = {} if plot_params is None else plot_params
    fig = plot_simulation_results(data, metrics, palette=palette, **plot_params)
    save_plot(fig, figure_name, output_dir)
    plt.close(fig)


def plot_figureS5(
    figure2_csv, euc_csv, corr_csv, metrics, palette, figure_name, output_dir
):
    """
    Combine CSVs for different distance metrics, plot simulation results, and save the figure.
    """
    distances_df = pd.read_csv(figure2_csv)
    euclidean_df = pd.read_csv(euc_csv)
    correlation_df = pd.read_csv(corr_csv)

    # Rename and add columns for the different distance metrics.
    distances_df.rename(columns={"mAP p-value": "mAP p-value (cosine)"}, inplace=True)
    distances_df["mAP p-value (Euclidean)"] = euclidean_df["mAP p-value"]
    distances_df["mAP p-value (correlation)"] = correlation_df["mAP p-value"]

    fig = plot_simulation_results(distances_df, metrics=metrics, palette=palette)
    if fig is None:
        fig = plt.gcf()
    save_plot(fig, figure_name, output_dir)
    plt.close(fig)


def main():
    output_dir = Path("outputs/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    figure2_palette = get_figure2_palette()
    figure_S5_palette = get_figureS5_palette()

    plot_params = None
    set_plotting_style(font_size=16, linewidth=1)
    sns.set_theme(style="white", rc={"lines.markeredgewidth": 0})

    # # uncomment for single-column format
    # set_plotting_style(font_size=5, linewidth=0.75)
    # plot_params = {
    #     "height": 1,
    #     "ylim": (0, 1),
    #     "font_size": 5,
    #     "markersize": 3,
    #     "labelrotation": 0,
    # }

    # Figure 2: Benchmarking mAP vs mp-value, MMD, and k-means on simulated data.
    plot_and_save(
        csv_path="results/fig2_results.csv",
        metrics=["mAP p-value", "mp-value", "MMD", "k-means"],
        palette=figure2_palette,
        figure_name="Figure2",
        output_dir=output_dir,
        plot_params=plot_params,
    )

    # Figure S2: Same as Figure 2, but with X axis spanning 0-100% feature change.
    plot_and_save(
        csv_path="results/fig_s2_results.csv",
        metrics=["mAP p-value", "mp-value", "MMD", "k-means"],
        palette=figure2_palette,
        figure_name="FigureS2",
        output_dir=output_dir,
    )

    # Figure S3: Benchmarking on simulated data with perturbed features ~ N(2,2).
    plot_and_save(
        csv_path="results/fig_s3_results.csv",
        metrics=["mAP p-value", "mp-value", "MMD", "k-means"],
        palette=figure2_palette,
        figure_name="FigureS3",
        output_dir=output_dir,
    )

    # Figure S4: Benchmarking on simulated data with Cauchy distributions.
    plot_and_save(
        csv_path="results/fig_s4_results.csv",
        metrics=["mAP p-value", "mp-value", "MMD", "k-means"],
        palette=figure2_palette,
        figure_name="FigureS4",
        output_dir=output_dir,
    )

    # Figure S5: Benchmarking mAP using different underlying distances.
    plot_figureS5(
        figure2_csv="results/fig2_results.csv",
        euc_csv="results/fig_s5_euc_results.csv",
        corr_csv="results/fig_s5_corr_results.csv",
        metrics=[
            "mAP p-value (cosine)",
            "mAP p-value (Euclidean)",
            "mAP p-value (correlation)",
        ],
        palette=figure_S5_palette,
        figure_name="FigureS5",
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()
