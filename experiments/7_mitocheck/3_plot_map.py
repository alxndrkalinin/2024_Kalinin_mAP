#!/usr/bin/env python3
# coding: utf-8

# 3. Plot mAP for morphological class and gene retrieval using CellProfiler and DeepProfiler features

import warnings
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from map_utils.plot import plot_map_scatter_kde, set_plotting_style, save_plot
from mitocheck_plot_utils import get_umap_palette, plot_gene_vs_ms, plot_mc_cp_vs_dp

warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")


def load_and_prepare_ap_results(ap_csv_path: str) -> pd.DataFrame:
    """
    Load AP results from CSV, rename the AP column, and update the Features values.
    """
    df = pd.read_csv(ap_csv_path)
    df = df.query("Features != 'cp_dp'")
    df.rename(columns={"average_precision": "AP"}, inplace=True)
    df["Features"] = df["Features"].replace(
        {"cp": "CellProfiler features", "dp": "DeepProfiler features"}
    )
    return df


def plot_map_results(ap_df: pd.DataFrame, fig_name: str, fig_dir: Path) -> None:
    """
    Plot the map results using the provided AP DataFrame and save the figure.
    """
    fig = plot_map_scatter_kde(
        ap_df,
        "Features",
        "",
        metric="AP",
        row=None,
        legend_loc="lower right",
        aspect=None,
        adjust=None,
        pr_x=0.5,
        pr_y=0.02,
        l_x=1.15,
        l_y=0.5,
        m_x=0.52,
        m_y=0.01,
    )
    save_plot(fig, fig_name, fig_dir)
    plt.close(fig)


def prepare_and_plot_gene_vs_ms(
    gene_ap: pd.DataFrame,
    morphoclass_ap: pd.DataFrame,
    fig_name_prefix: str,
    fig_dir: Path,
    plot_kwargs: dict = None,
) -> None:
    """
    Prepare DataFrames for gene vs. morphological class retrieval and plot using plot_gene_vs_ms.
    Two plots are generated—one for CellProfiler and one for DeepProfiler features.
    """
    # For CellProfiler features.
    aps_cp = pd.DataFrame(
        {
            "AP, Gene": gene_ap.query("Features == 'CellProfiler features'")["AP"],
            "Gene p < 0.05": gene_ap.query("Features == 'CellProfiler features'")[
                "p < 0.05"
            ],
            "AP, Morphological Class (MC)": morphoclass_ap.query(
                "Features == 'CellProfiler features'"
            )["AP"],
            "MC p < 0.05": morphoclass_ap.query("Features == 'CellProfiler features'")[
                "p < 0.05"
            ],
        }
    )
    fig_cp = plot_gene_vs_ms(aps_cp, title="CellProfiler features", **plot_kwargs)
    save_plot(fig_cp, f"{fig_name_prefix}_CellProfiler", fig_dir)

    # For DeepProfiler features.
    aps_dp = pd.DataFrame(
        {
            "AP, Gene": gene_ap.query("Features == 'DeepProfiler features'")["AP"],
            "Gene p < 0.05": gene_ap.query("Features == 'DeepProfiler features'")[
                "p < 0.05"
            ],
            "AP, Morphological Class (MC)": morphoclass_ap.query(
                "Features == 'DeepProfiler features'"
            )["AP"],
            "MC p < 0.05": morphoclass_ap.query("Features == 'DeepProfiler features'")[
                "p < 0.05"
            ],
        }
    )
    fig_dp = plot_gene_vs_ms(aps_dp, title="DeepProfiler features", **plot_kwargs)
    save_plot(fig_dp, f"{fig_name_prefix}_DeepProfiler", fig_dir)


def prepare_and_plot_mc_cp_vs_dp(
    morphoclass_ap: pd.DataFrame,
    fig_name: str,
    fig_dir: Path,
    plot_kwargs: dict = None,
) -> None:
    """
    Prepare a DataFrame comparing CP and DP retrieval for each morphological class and plot using plot_mc_cp_vs_dp.
    """
    aps = pd.DataFrame(
        {
            "Morphological Class": morphoclass_ap.query(
                "Features == 'CellProfiler features'"
            )["Mitocheck_Phenotypic_Class"].reset_index(drop=True),
            "AP, CellProfiler features": morphoclass_ap.query(
                "Features == 'CellProfiler features'"
            )["AP"].reset_index(drop=True),
            "CellProfiler retrieved": morphoclass_ap.query(
                "Features == 'CellProfiler features'"
            )["p < 0.05"].reset_index(drop=True),
            "AP, DeepProfiler features": morphoclass_ap.query(
                "Features == 'DeepProfiler features'"
            )["AP"].reset_index(drop=True),
            "DeepProfiler retrieved": morphoclass_ap.query(
                "Features == 'DeepProfiler features'"
            )["p < 0.05"].reset_index(drop=True),
        }
    )
    plot_mc_cp_vs_dp(aps, hue="Morphological Class", **plot_kwargs)
    fig = plt.gcf()
    save_plot(fig, fig_name, fig_dir)
    plt.close(fig)


def plot_boxplot_comparison(
    morphoclass_ap: pd.DataFrame, fig_name: str, fig_dir: Path
) -> None:
    """
    Create a boxplot comparing AP (for profiles with p < 0.05) across morphological classes.
    """
    cp_means = (
        morphoclass_ap.query("Features == 'CellProfiler features' and `p < 0.05`")
        .groupby("Mitocheck_Phenotypic_Class")["AP"]
        .mean()
    )
    dp_means = (
        morphoclass_ap.query("Features == 'DeepProfiler features' and `p < 0.05`")
        .groupby("Mitocheck_Phenotypic_Class")["AP"]
        .mean()
    )
    sorted_order = (cp_means - dp_means).sort_values(ascending=False).index

    df_plot = morphoclass_ap.copy()
    df_plot["Mitocheck_Phenotypic_Class"] = pd.Categorical(
        df_plot["Mitocheck_Phenotypic_Class"], categories=sorted_order, ordered=True
    )
    df_plot = df_plot.sort_values("Mitocheck_Phenotypic_Class")
    fig = plt.figure(figsize=(12, 6))
    ax = sns.boxplot(
        data=df_plot[df_plot["p < 0.05"]].rename(
            columns={
                "AP": "AP, p < 0.05",
                "Mitocheck_Phenotypic_Class": "Morphological Class",
            }
        ),
        x="Morphological Class",
        y="AP, p < 0.05",
        hue="Features",
        fill=False,
        palette="Set2",
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.legend(frameon=False, title=None, loc="upper center")
    save_plot(fig, fig_name, fig_dir)
    plt.close(fig)


def prepare_aps_retrieval_classes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare a DataFrame with retrieval info (Cell ID, morphological class, AP and retrieval flags)
    for both CellProfiler and DeepProfiler features.
    """
    aps_rc = pd.DataFrame(
        {
            "Cell ID": df.query("Features == 'CellProfiler features'")[
                "Cell_UUID"
            ].reset_index(drop=True),
            "Morphological Class": df.query("Features == 'CellProfiler features'")[
                "Mitocheck_Phenotypic_Class"
            ].reset_index(drop=True),
            "AP, CellProfiler features": df.query(
                "Features == 'CellProfiler features'"
            )["AP"].reset_index(drop=True),
            "CellProfiler retrieved": df.query("Features == 'CellProfiler features'")[
                "p < 0.05"
            ].reset_index(drop=True),
            "AP, DeepProfiler features": df.query(
                "Features == 'DeepProfiler features'"
            )["AP"].reset_index(drop=True),
            "DeepProfiler retrieved": df.query("Features == 'DeepProfiler features'")[
                "p < 0.05"
            ].reset_index(drop=True),
        }
    )
    return aps_rc


def load_sc_profiles() -> pd.DataFrame:
    """
    Load the preprocessed, merged single-cell profiles from the merged file.
    This file already contains both training and negative control profiles.
    """
    df = pd.read_parquet("outputs/training_sc_fs.parquet")
    df = df.query(
        "Metadata_Gene != 'failed QC' and Mitocheck_Phenotypic_Class != 'OutOfFocus'"
    )
    return df


def assign_umap_hue_and_order(
    df: pd.DataFrame, prefix: str, retrieved_col: str
) -> pd.DataFrame:
    """
    Rename UMAP columns based on the given prefix (CP or DP),
    assign a 'hue' based on the retrieved flag, and sort by a plot_order column.
    """
    df = df.copy()
    if prefix == "CP":
        rename_map = {"CP_umap_0": "UMAP_1", "CP_umap_1": "UMAP_2"}
    else:
        rename_map = {"DP_umap_0": "UMAP_1", "DP_umap_1": "UMAP_2"}
    df.rename(columns=rename_map, inplace=True)

    df["hue"] = "Negative control"
    df.loc[df["class"] & ~df[retrieved_col], "hue"] = "False"
    df.loc[df["class"] & df[retrieved_col], "hue"] = "True"

    order = {"Negative control": 0, "False": 1, "True": 2}
    df["plot_order"] = df["hue"].map(order)
    return df.sort_values("plot_order")


def plot_umap_scatter(ax, data: pd.DataFrame, title: str) -> None:
    """
    Create a scatterplot on the given axis using UMAP data.
    """
    sns.scatterplot(
        data=data,
        x="UMAP_1",
        y="UMAP_2",
        hue="hue",
        hue_order=["True", "False", "Negative control"],
        palette=get_umap_palette(n=data.hue.nunique())[::-1],
        s=20,
        ax=ax,
    )
    ax.set_title(title)
    ax.get_legend().remove()


def plot_umap_by_morpho_class(aps_retrieval, fig_dir, output_dir):
    """
    Load UMAP embeddings from disk, merge with retrieval info, and create scatter plots.
    Saves each figure.
    """
    for morpho_class in aps_retrieval["Morphological Class"].unique():
        umap_path = output_dir / f"umap_{morpho_class}.parquet"
        if not umap_path.exists():
            print(f"Warning: {umap_path} not found. Skipping {morpho_class}.")
            continue

        embed_df = pd.read_parquet(umap_path)

        merged_data = embed_df.merge(
            aps_retrieval, how="left", left_on=["Cell_UUID"], right_on=["Cell ID"]
        )

        cp_retrieval = merged_data["CellProfiler retrieved"].mean()
        dp_retrieval = merged_data["DeepProfiler retrieved"].mean()

        merged_data.fillna(
            {"CellProfiler retrieved": False, "DeepProfiler retrieved": False},
            inplace=True,
        )
        merged_data["class"] = merged_data["Mitocheck_Phenotypic_Class"] == morpho_class

        merged_data_cp = assign_umap_hue_and_order(
            merged_data, "CP", "CellProfiler retrieved"
        )
        merged_data_dp = assign_umap_hue_and_order(
            merged_data, "DP", "DeepProfiler retrieved"
        )

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))
        plot_umap_scatter(
            ax[0], merged_data_cp, f"CellProfiler features ({cp_retrieval:.0%})"
        )
        plot_umap_scatter(
            ax[1], merged_data_dp, f"DeepProfiler features ({dp_retrieval:.0%})"
        )

        fig.suptitle(f"Morphological Class = {morpho_class}")
        handles, labels = ax[1].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            title="p < 0.05",
            loc="upper center",
            bbox_to_anchor=(1.0, 0.5),
            frameon=False,
        )

        save_plot(fig, f"UMAP_{morpho_class}", fig_dir)
        plt.close(fig)


def main():
    output_dir = Path("outputs")
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    gene_ap = load_and_prepare_ap_results(output_dir / "Metadata_Gene_ap_results.csv")
    morphoclass_ap = load_and_prepare_ap_results(
        output_dir / "Mitocheck_Phenotypic_Class_ap_results.csv"
    )

    # Plot map results for both gene and morphological class AP.
    set_plotting_style(font_size=16, linewidth=1)
    plot_map_results(gene_ap, "Map_Gene_AP", fig_dir)
    plot_map_results(morphoclass_ap, "Map_Morphoclass_AP", fig_dir)

    # Plot gene vs. morphological class AP retrieval.
    set_plotting_style(font_size=5, linewidth=0.35)
    prepare_and_plot_gene_vs_ms(
        gene_ap,
        morphoclass_ap,
        "Gene_vs_Morphoclass",
        fig_dir,
        plot_kwargs={
            "width": 1.2,
            "height": 1.2,
            "point_size": 1,
            "title_font_size": 5,
            "legend_markersize": 2,
        },
    )

    # Plot morphological class comparison between CP and DP features.
    prepare_and_plot_mc_cp_vs_dp(
        morphoclass_ap,
        "MC_CP_vs_DP",
        fig_dir,
        plot_kwargs={
            "width": 1.3,
            "height": 1.3,
            "point_size": 1.5,
            "legend_markersize": 2,
        },
    )

    # Plot a boxplot comparison by morphological class.
    set_plotting_style(font_size=16, linewidth=1)
    plot_boxplot_comparison(morphoclass_ap, "Boxplot_Morphoclass", fig_dir)

    # Plot UMAP embeddings by morphological class.
    aps_retrieval = prepare_aps_retrieval_classes(morphoclass_ap)
    plot_umap_by_morpho_class(aps_retrieval, fig_dir, output_dir)


if __name__ == "__main__":
    main()
