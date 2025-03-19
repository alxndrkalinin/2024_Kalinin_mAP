from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from map_utils.plot import plot_map_scatter, set_plotting_style


def save_figure(fig, filename_prefix: str, output_dir: Path) -> None:
    """
    Save a matplotlib figure in both PNG and SVG formats with tight layout.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    png_file = output_dir / f"{filename_prefix}.png"
    svg_file = output_dir / f"{filename_prefix}.svg"
    fig.savefig(png_file, format="png", bbox_inches="tight")
    fig.savefig(svg_file, format="svg", bbox_inches="tight")
    print(f"Saved figure as:\n  {png_file}\n  {svg_file}")


def process_and_plot_map(
    input_file: str, const_value: str, pr_x: float, pr_y: float, l_x: float, l_y: float
):
    """
    Load a map CSV file, compute -log10(mAP p-value), rename columns,
    add a constant column, print a preview, and generate a plot using plot_map_scatter.

    Returns:
        A matplotlib Figure.
    """
    df = pd.read_csv(input_file)
    df["-log10(mAP p-value)"] = -np.log10(df.corrected_p_value)
    df.rename(
        {"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        axis=1,
        inplace=True,
    )
    print(f"Preview of {input_file}:")
    print(df.head(10))
    df["const"] = const_value
    fig = plot_map_scatter(df, "const", "", pr_x=pr_x, pr_y=pr_y, l_x=l_x, l_y=l_y)
    if fig is None:
        fig = plt.gcf()
    return fig


def main():
    set_plotting_style()
    output_dir = Path("outputs")

    # Process and plot technical map
    technical_file = "outputs/map_technical.csv.gz"
    print("Processing technical map...")
    fig_tech = process_and_plot_map(
        technical_file,
        const_value="dataset: cpg0004",
        pr_x=0.45,
        pr_y=0.02,
        l_x=1.05,
        l_y=0.58,
    )
    save_figure(fig_tech, "technical_map", output_dir)

    # Calculate biological mAP using lincs.mean_average_precision
    from lincs import mean_average_precision

    print("Calculating biological mAP...")
    mean_average_precision(
        "outputs/ap_biological.csv.gz",
        "outputs/map_biological.csv.gz",
        sameby=["Metadata_target"],
        pvalue_threshold=0.05,
    )

    # Process and plot biological map
    biological_file = "outputs/map_biological.csv.gz"
    print("Processing biological map...")
    df_bio = pd.read_csv(biological_file)
    df_bio["-log10(mAP p-value)"] = -np.log10(df_bio.corrected_p_value)
    df_bio.rename(
        {"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        axis=1,
        inplace=True,
    )
    print("Preview of biological map:")
    print(df_bio.head(10))
    unique_targets = df_bio.Metadata_target.nunique()
    print(f"Number of unique Metadata_target entries: {unique_targets}")
    df_bio["const"] = "dataset: cpg0004"
    fig_bio = plot_map_scatter(
        df_bio, "const", "", pr_x=0.45, pr_y=0.02, l_x=1.1, l_y=0.58
    )
    if fig_bio is None:
        fig_bio = plt.gcf()
    save_figure(fig_bio, "biological_map", output_dir)


if __name__ == "__main__":
    main()
