from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from map_utils.plot import (
    add_corner_text_annotations,
    plot_map_scatter,
    set_plotting_style,
)


def load_and_combine_activity_maps() -> pd.DataFrame:
    """
    Load and combine phenotypic activity maps for Cell Painting and nELISA.
    """
    cp_act = pd.read_csv("outputs/cp_map_activity_results.csv")
    ne_act = pd.read_csv("outputs/ne_map_activity_results.csv")
    print(
        "Cell Painting activity shape:",
        cp_act.shape,
        "nELISA activity shape:",
        ne_act.shape,
    )

    cp_act["Assay"] = "Cell Painting"
    ne_act["Assay"] = "nELISA"

    tech_map = pd.concat([cp_act, ne_act], ignore_index=True)
    tech_map["-log10(mAP p-value)"] = -np.log10(tech_map["corrected_p_value"])
    tech_map.rename(
        columns={"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        inplace=True,
    )
    print("Combined activity map shape:", tech_map.shape)
    return tech_map


def load_and_combine_consistency_maps() -> pd.DataFrame:
    """
    Load and combine phenotypic consistency maps for Cell Painting and nELISA.
    """
    cp_cons = pd.read_csv("outputs/cp_all_map_consistency_results.csv")
    ne_cons = pd.read_csv("outputs/ne_all_map_consistency_results.csv")
    print(
        "Cell Painting consistency shape:",
        cp_cons.shape,
        "nELISA consistency shape:",
        ne_cons.shape,
    )

    cp_cons["Assay"] = "Cell Painting"
    ne_cons["Assay"] = "nELISA"

    bio_map_all = pd.concat([cp_cons, ne_cons], ignore_index=True)
    bio_map_all["-log10(mAP p-value)"] = -np.log10(bio_map_all["corrected_p_value"])
    bio_map_all.rename(
        columns={"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        inplace=True,
    )
    print("Combined consistency map shape:", bio_map_all.shape)
    return bio_map_all


def plot_scatter_comparison(
    data: pd.DataFrame,
    pivot_index: str,
    x_label: str,
    y_label: str,
    threshold: float,
    threshold_text_v_offset: float,
    corner_h_offset: float,
    corner_v_offset: float,
    save_prefix: str = None,
    output_dir: Path = None,
) -> None:
    """
    Create a scatter plot comparing -log10(mAP p-value) for two assays and optionally save the figure.

    Parameters
    ----------
    data : pd.DataFrame
        Combined DataFrame with a column "Assay" (values: "Cell Painting" and "nELISA")
        and a "-log10(mAP p-value)" column.
    pivot_index : str
        Column name to pivot on (e.g., "Metadata_broad_sample" or "Metadata_target").
    x_label : str
        Label for the x-axis (assumed to be the value for Cell Painting).
    y_label : str
        Label for the y-axis (assumed to be the value for nELISA).
    threshold : float
        Threshold value to draw horizontal and vertical lines.
    threshold_text_v_offset : float
        Vertical offset for the threshold text.
    corner_h_offset : float
        Horizontal offset for the corner text annotations.
    corner_v_offset : float
        Vertical offset for the corner text annotations.
    save_prefix : str, optional
        If provided, the figure will be saved with this prefix.
    output_dir : Path, optional
        Directory where the figure will be saved (must be provided if save_prefix is given).
    """
    df_pivot = data.pivot(
        index=pivot_index, columns="Assay", values="-log10(mAP p-value)"
    ).reset_index()

    fig = plt.figure(figsize=(4.5, 4.5))
    ax = sns.scatterplot(data=df_pivot, x=x_label, y=y_label)

    ax.axhline(threshold, color="grey", linestyle="--")
    ax.axvline(threshold, ymax=0.9, color="grey", linestyle="--")
    ax.text(
        0.02,
        threshold / df_pivot[y_label].max() + threshold_text_v_offset,
        "p=0.05",
        transform=ax.transAxes,
        color="grey",
        fontsize=12,
        fontstyle="italic",
    )
    plt.xlabel(f"-log10(mAP p-value), {x_label}")
    plt.ylabel(f"-log10(mAP p-value), {y_label}")

    lim = max(df_pivot[x_label].max(), df_pivot[y_label].max())
    plt.plot([0, lim], [0, lim], "r--")

    add_corner_text_annotations(
        ax,
        df_pivot,
        "nELISA",
        "Cell Painting",
        prefix="Retrieved: ",
        h_offset=corner_h_offset,
        v_offset=corner_v_offset,
    )

    if save_prefix is not None and output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        png_file = output_dir / f"{save_prefix}.png"
        svg_file = output_dir / f"{save_prefix}.svg"
        fig.savefig(png_file, format="png", bbox_inches="tight")
        fig.savefig(svg_file, format="svg", bbox_inches="tight")
        print(f"Saved scatter plot as:\n  {png_file}\n  {svg_file}")

    plt.show()


def save_map_plot(fig, filename_prefix: str, output_dir: Path) -> None:
    """
    Save a given figure in both PNG and SVG formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    filename_prefix : str
        Prefix for the saved files.
    output_dir : Path
        Directory where files will be saved.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    png_file = output_dir / f"{filename_prefix}.png"
    svg_file = output_dir / f"{filename_prefix}.svg"
    fig.savefig(png_file, format="png", bbox_inches="tight")
    fig.savefig(svg_file, format="svg", bbox_inches="tight")
    print(f"Saved map plot as:\n  {png_file}\n  {svg_file}")
    plt.close(fig)


def main():
    set_plotting_style()
    output_dir = Path("outputs")

    # Load and plot phenotypic activity maps
    tech_map = load_and_combine_activity_maps()
    print("Plotting phenotypic activity map...")
    fig1 = plot_map_scatter(
        tech_map, "Assay", "", pr_x=0.5, pr_y=0.02, m_x=0.52, m_y=0.02, kde_y=0.75
    )
    if fig1 is None:
        fig1 = plt.gcf()
    save_map_plot(fig1, "activity_map", output_dir)

    plot_scatter_comparison(
        data=tech_map,
        pivot_index="Metadata_broad_sample",
        x_label="Cell Painting",
        y_label="nELISA",
        threshold=-np.log10(0.05),
        threshold_text_v_offset=0.02,
        corner_h_offset=0.02,
        corner_v_offset=0.01,
        save_prefix="activity_scatter",
        output_dir=output_dir,
    )

    # Load and plot phenotypic consistency maps
    bio_map_all = load_and_combine_consistency_maps()
    print("Plotting phenotypic consistency map...")
    fig2 = plot_map_scatter(
        bio_map_all, "Assay", "", pr_x=0.5, pr_y=0.02, m_x=0.52, m_y=0.02, kde_y=0.65
    )
    if fig2 is None:
        fig2 = plt.gcf()
    save_map_plot(fig2, "consistency_map", output_dir)

    plot_scatter_comparison(
        data=bio_map_all,
        pivot_index="Metadata_target",
        x_label="Cell Painting",
        y_label="nELISA",
        threshold=-np.log10(0.05),
        threshold_text_v_offset=-0.075,
        corner_h_offset=0.06,
        corner_v_offset=0.02,
        save_prefix="consistency_scatter",
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()
