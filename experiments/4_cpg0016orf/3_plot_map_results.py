#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from map_utils.plot import plot_map_scatter_kde, set_plotting_style


def process_map(file_path: Path) -> pd.DataFrame:
    """
    Load and process a map file (Parquet) and add computed columns.
    """
    df = pd.read_parquet(file_path)
    df["-log10(mAP p-value)"] = -np.log10(df["corrected_p_value"])
    df.rename(
        {"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        axis=1,
        inplace=True,
    )
    df["const"] = "dataset: cpg0016[orf]"
    print(f"Processed map shape for {file_path}: {df.shape}")
    return df


def plot_and_save_map(data: pd.DataFrame, prefix: str, output_dir: Path) -> None:
    """
    Plot the map using plot_map_scatter and save the figure in PNG and SVG formats.
    """
    fig = plot_map_scatter_kde(
        data, "const", "", pr_x=0.45, pr_y=0.02, l_x=1.1, l_y=0.58
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    png_file = output_dir / f"{prefix}.png"
    svg_file = output_dir / f"{prefix}.svg"

    fig.savefig(png_file, format="png", bbox_inches="tight")
    fig.savefig(svg_file, format="svg", bbox_inches="tight")
    print(f"Saved {prefix} plot as:\n  {png_file}\n  {svg_file}")
    plt.close(fig)


def main():
    set_plotting_style()

    output_dir = Path("outputs")
    activity_file = (
        output_dir
        / "wellpos_cc_var_mad_outlier_featselect_sphering_harmony_activity.parquet"
    )
    consistency_file = output_dir / "corum_complex_consistency.parquet"

    activity_map = process_map(activity_file)
    plot_and_save_map(activity_map, "activity_map", output_dir)

    consistency_map = process_map(consistency_file)
    plot_and_save_map(consistency_map, "consistency_map", output_dir)


if __name__ == "__main__":
    main()
