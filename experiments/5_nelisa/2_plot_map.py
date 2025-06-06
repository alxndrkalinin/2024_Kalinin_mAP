from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from map_utils.plot import (
    save_plot,
    plot_map_scatter_kde,
    set_plotting_style,
    add_corner_text_annotations,
)


def load_and_combine_maps(cp_path: str, ne_path: str, label: str) -> pd.DataFrame:
    """
    Load and combine Cell Painting and nELISA results.

    Parameters:
    - cp_path: Path to Cell Painting CSV file.
    - ne_path: Path to nELISA CSV file.
    - label: Descriptive label for printing (e.g., "activity", "consistency").

    Returns:
    - Combined DataFrame with assay labels and standardized columns.
    """
    cp = pd.read_csv(cp_path)
    ne = pd.read_csv(ne_path)

    print(f"Cell Painting {label} shape:", cp.shape, f"nELISA {label} shape:", ne.shape)

    cp["Assay"] = "Cell Painting"
    ne["Assay"] = "nELISA"

    df = pd.concat([cp, ne], ignore_index=True)
    df["-log10(mAP p-value)"] = -np.log10(df["corrected_p_value"])
    df.rename(
        columns={"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        inplace=True,
    )
    print(f"Combined {label} map shape:", df.shape)
    return df


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
    point_size: int = 5,
    h_offset_bottom: float = None,
    v_offset_bottom: float = None,
    fig_size: tuple = (4.5, 4.5),
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

    fig = plt.figure(figsize=fig_size)
    ax = sns.scatterplot(data=df_pivot, x=x_label, y=y_label, s=point_size)

    ax.axhline(threshold, color="grey", linestyle="--")
    ax.axvline(threshold, ymax=0.9, color="grey", linestyle="--")
    # ax.text(
    #     0.02,
    #     threshold / df_pivot[y_label].max() + threshold_text_v_offset,
    #     "p=0.05",
    #     transform=ax.transAxes,
    #     color="grey",
    #     fontsize=12,
    #     fontstyle="italic",
    # )
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
        h_offset_bottom=h_offset_bottom,
        v_offset_bottom=v_offset_bottom,
    )
    save_plot(fig, save_prefix, output_dir)
    plt.show()


def main():
    set_plotting_style(font_size=5, linewidth=0.35)
    output_dir = Path("outputs")

    # Load and plot phenotypic activity maps
    activity_df = load_and_combine_maps(
        "outputs/cp_map_activity_results.csv",
        "outputs/ne_map_activity_results.csv",
        label="activity",
    )

    print("Plotting phenotypic activity map...")
    fig1 = plot_map_scatter_kde(
        activity_df,
        "Assay",
        "",
        pr_x=0.45,
        pr_y=0.02,
        m_x=0.52,
        m_y=0.02,
        l_x=1.1,
        l_y=0.575,
        kde_y=0.75,
        size_rescale=0.31,
        point_size=5,
        legend=True,
        legend_frameon=True,
    )
    save_plot(fig1, "activity_map", output_dir)

    plot_scatter_comparison(
        data=activity_df,
        pivot_index="Metadata_broad_sample",
        x_label="Cell Painting",
        y_label="nELISA",
        threshold=-np.log10(0.05),
        threshold_text_v_offset=0.02,
        corner_h_offset=0.02,
        corner_v_offset=0.01,
        save_prefix="activity_scatter",
        output_dir=output_dir,
        fig_size=(1, 1),
        point_size=3,
        h_offset_bottom=0,
        v_offset_bottom=-0.01,
    )

    # Load and plot phenotypic consistency maps
    consistency_df = load_and_combine_maps(
        "outputs/cp_all_map_consistency_results.csv",
        "outputs/ne_all_map_consistency_results.csv",
        label="consistency",
    )
    print("Plotting phenotypic consistency map...")
    fig2 = plot_map_scatter_kde(
        consistency_df,
        "Assay",
        "",
        pr_x=0.5,
        pr_y=0.02,
        m_x=0.52,
        m_y=0.02,
        kde_y=0.65,
        size_rescale=0.31,
        point_size=5,
    )
    save_plot(fig2, "consistency_map", output_dir)

    plot_scatter_comparison(
        data=consistency_df,
        pivot_index="Metadata_target",
        x_label="Cell Painting",
        y_label="nELISA",
        threshold=-np.log10(0.05),
        threshold_text_v_offset=-0.075,
        corner_h_offset=0.06,
        corner_v_offset=0.02,
        save_prefix="consistency_scatter",
        output_dir=output_dir,
        fig_size=(1, 1),
        point_size=3,
        h_offset_bottom=0.1,
        v_offset_bottom=-0.01,
    )


if __name__ == "__main__":
    main()
