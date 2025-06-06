# 4. Plot mAP vs relative activity for bulk and single-cell data (Figs 4C & 5A)
#
# The perturbseq data come from the CRISPRi experiment from:
#
# > Jost, M., Santos, D.A., Saunders, R.A. et al. Titrating gene expression using libraries of systematically attenuated CRISPR guide RNAs. Nat Biotechnol 38, 355–364 (2020). https://doi.org/10.1038/s41587-019-0387-5
#
# and relative activity is defined as:
#
# The fold-knockdown of each mismatched variant divided by the fold-knockdown of the perfectly-matched sgRNA.

from pathlib import Path

import numpy as np
import pandas as pd
import plotnine as gg
import seaborn as sns

from map_utils.plot import set_plotting_style, save_plot, plot_map_activity_scatter


def load_bulk_map(gse_id: str, results_dir: Path) -> pd.DataFrame:
    """
    Load bulk mAP results from file and process the DataFrame.
    """
    bulk_file = results_dir / f"{gse_id}_map.tsv"
    bulk_df = pd.read_csv(bulk_file, sep="\t")
    bulk_df["gene"] = pd.Categorical(bulk_df.gene, categories=bulk_df.gene.unique())
    print("Bulk map data shape:", bulk_df.shape)
    print(bulk_df.head(2))
    return bulk_df


def make_colormap() -> dict:
    """
    Create a color map for the scatter plot.
    """
    palette = sns.color_palette().as_hex()
    return {False: palette[0], True: palette[1]}


def create_bulk_scatter_plots(
    bulk_df: pd.DataFrame,
    output_dir: Path,
    gse_id: str,
    colormap=None,
    scatterplot_size=(3, 2),
    faceted_size=(5, 6),
):
    """
    Create and save global and facet scatter plots for bulk data.
    """
    global_gg = plot_map_activity_scatter(
        df=bulk_df,
        x="relative_activity_day5",
        y="mAP",
        aes_mapping={"fill": "p < 0.05"},
        add_linear_model=True,
        x_label="Mismatched guide RNA activity relative to perfect match",
        y_label="mAP",
        scale_manual_values=colormap,
        point_size=1.0,
        point_alpha=1.0,
        y_range=(0, 1.0),
        x_range=(0, bulk_df.relative_activity_day5.max()),
        x_scale={
            "breaks": [0, 0.5, 1],
            "labels": ["0", "0.5", "1"],
            "minor_breaks": [],
        },
        y_scale={
            "breaks": [0, 0.5, 1],
            "labels": ["0", "0.5", "1"],
            "minor_breaks": [],
        },
        width=scatterplot_size[0],
        height=scatterplot_size[1],
    )

    save_plot(
        global_gg,
        f"{gse_id}_crispri_map_relative_activity_global_bulk",
        output_dir,
        dpi=500,
        width=scatterplot_size[0],
        height=scatterplot_size[1],
    )

    gene_subset = bulk_df.gene.unique()[:]
    gene_gg = plot_map_activity_scatter(
        df=bulk_df.query("gene in @gene_subset"),
        x="relative_activity_day5",
        y="mAP",
        facet="gene",
        facet_ncol=5,
        aes_mapping={"fill": "p < 0.05"},
        x_label="Mismatched guide RNA activity relative to perfect match",
        y_label="mAP",
        y_range=(0, 1.05),
        x_range=(0, bulk_df.relative_activity_day5.max()),
        scale_manual_values=colormap,
        point_size=1.0,
        point_alpha=0.8,
        point_stroke=0.05,
        guide_title="mAP",
        x_scale={
            "breaks": [0, 0.5, 1],
            "labels": ["0", "0.5", "1"],
            "minor_breaks": [],
        },
        y_scale={"breaks": [0, 1], "labels": ["0", "1"], "minor_breaks": []},
        width=faceted_size[0],
        height=faceted_size[1],
    )
    save_plot(
        gene_gg,
        f"{gse_id}_crispri_map_relative_activity_facet_bulk",
        output_dir,
        dpi=500,
    )


def load_single_cell_data(
    gse_id: str, data_dir: Path, results_dir: Path, bulk_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Load finalized single cell data, mAP results, and UMAP embeddings, then merge them with bulk data.
    """
    # Load finalized single cell gene expression data
    gene_exp_file = results_dir / f"{gse_id}_final_analytical.tsv.gz"
    sc_gene_exp_df = pd.read_csv(gene_exp_file, sep="\t", low_memory=False)
    # Get gene features (columns not starting with "Metadata_")
    gene_features = [
        col for col in sc_gene_exp_df.columns if not col.startswith("Metadata_")
    ]
    print("Single-cell gene expression data shape:", sc_gene_exp_df.shape)

    # Filter single cell data
    sc_gene_exp_df = (
        sc_gene_exp_df.query("Metadata_gene_identity != '*'")
        .query("Metadata_good_coverage == True")
        .reset_index(drop=True)
    )
    sc_gene_exp_df["Metadata_reference_index"] = np.where(
        sc_gene_exp_df["Metadata_gene_identity"] == "neg", sc_gene_exp_df.index, -1
    )

    # Load single cell mAP results
    sc_map_file = results_dir / f"{gse_id}_single_cell_map.tsv.gz"
    sc_map_df = pd.read_csv(sc_map_file, sep="\t")
    sc_map_df["-log10(mAP p-value)"] = -np.log10(sc_map_df["p_value"])
    sc_map_df.rename(columns={"mean_ap": "mAP"}, inplace=True)
    sc_map_df["p < 0.05"] = sc_map_df["p_value"] < 0.05

    # Load UMAP embeddings
    sc_umap_file = results_dir / f"{gse_id}_single_cell_umap_embeddings.tsv.gz"
    sc_umap_df = pd.read_csv(sc_umap_file, sep="\t")

    # Merge UMAP embeddings with mAP results
    sc_df = sc_umap_df.merge(
        sc_map_df,
        on=[
            "Metadata_cell_identity",
            "Metadata_guide_identity",
            "Metadata_gene_identity",
        ],
        how="right",
    )
    # Merge with bulk data to bring in the 'gene' column
    sc_df = sc_df.merge(
        bulk_df,
        on=["Metadata_guide_identity", "Metadata_gene_identity"],
        how="outer",
        suffixes=["", "_bulk_activity"],
    )
    sc_df["gene"] = pd.Categorical(sc_df.gene, categories=bulk_df.gene.unique())
    sc_df["Metadata_gene_identity"] = pd.Categorical(
        sc_df.Metadata_gene_identity,
        categories=["neg"] + bulk_df.gene.unique().tolist(),
    )
    print("Merged single-cell data shape:", sc_df.shape)
    return sc_df, gene_features


def create_singlecell_global_plot(
    sc_df: pd.DataFrame, output_dir: Path, gse_id: str, colormap=None, fig_size=(3, 2)
):
    """
    Create and save a global scatter plot (density and points) for single cell data.
    """
    singlecell_scatter = plot_map_activity_scatter(
        df=sc_df,
        x="relative_activity_day5",
        y="mAP",
        aes_mapping={"fill": "p < 0.05"},
        density_aes={"color": "p < 0.05"},
        density_size=0.1,
        x_label="Mismatched guide RNA activity relative to perfect match",
        y_label="mAP",
        scale_manual_values=colormap,
        point_size=0.3,
        point_alpha=0.3,
        guide_title="mAP",
        x_scale={
            "breaks": [0, 0.5, 1],
            "labels": ["0", "0.5", "1"],
            "minor_breaks": [],
        },
        y_scale={
            "breaks": [0, 0.5, 1],
            "labels": ["0", "0.5", "1"],
            "minor_breaks": [],
        },
        width=fig_size[0],
        height=fig_size[1],
    )
    save_plot(
        singlecell_scatter,
        f"{gse_id}_singlecell_global",
        output_dir,
        dpi=500,
        width=fig_size[0],
        height=fig_size[1],
    )


def create_singlecell_facet_plot(
    sc_df: pd.DataFrame, output_dir: Path, gse_id: str, colormap=None, fig_size=(3, 2)
):
    """
    Create and save a faceted scatter plot for single cell data by gene.
    """
    facet_gg = plot_map_activity_scatter(
        sc_df.dropna(subset=["gene"]),
        x="relative_activity_day5",
        y="mAP",
        facet="gene",
        facet_ncol=5,
        aes_mapping={"fill": "p < 0.05"},
        density_aes={"color": "p < 0.05"},
        density_size=0.1,
        point_size=0.1,
        point_alpha=0.05,
        x_label="Mismatched guide RNA activity relative to perfect match",
        y_label="mAP",
        y_range=(0, 1.05),
        x_range=(0, sc_df.relative_activity_day5.max()),
        scale_manual_values=colormap,
        x_scale={
            "breaks": [0, 0.5, 1],
            "labels": ["0", "0.5", "1"],
            "minor_breaks": [],
        },
        y_scale={"breaks": [0, 1], "labels": ["0", "1"], "minor_breaks": []},
        width=fig_size[0],
        height=fig_size[1],
    )
    save_plot(
        facet_gg, f"{gse_id}_singlecell_facet", output_dir, dpi=500, height=5, width=6
    )


def plot_individual_umaps(sc_df: pd.DataFrame, output_dir: Path, gse_id: str):
    """
    For each gene (excluding negatives and invalids), create a UMAP facet plot and save.
    """
    gene_umap_dir = output_dir / "gene_umaps"
    gene_umap_dir.mkdir(parents=True, exist_ok=True)

    valid_genes = [
        gene
        for gene in sc_df.gene.unique()
        if gene not in ["neg", "*", "nan"] and pd.notna(gene)
    ]

    for gene in valid_genes:
        print(f"Processing UMAP for gene: {gene}")
        gene_embedding_df = sc_df.query("Metadata_gene_identity == @gene")
        # First create a consistent string representation of relative_activity_day5
        rounded_values = gene_embedding_df.relative_activity_day5.round(3)
        activity_strings = rounded_values.apply(
            lambda x: f"{x:.3f}" if pd.notna(x) else "NaN"
        )

        # Create facet labels with consistent formatting
        gene_embedding_df = gene_embedding_df.assign(
            map_facet_label=gene_embedding_df.Metadata_gene_identity.astype(str)
            + " "
            + activity_strings
        )

        # Replace label for negatives
        gene_embedding_df.loc[
            gene_embedding_df.Metadata_gene_identity == "neg", "map_facet_label"
        ] = "Negative Ctrl"

        # Create facet_order with the same formatting
        unique_activities = sorted(rounded_values.dropna().unique())
        facet_order = ["Negative Ctrl"] + [
            f"{gene} {val:.3f}" for gene in [gene] for val in unique_activities
        ]

        gene_embedding_df["map_facet_label"] = pd.Categorical(
            gene_embedding_df.map_facet_label, categories=facet_order
        )

        gene_umap_gg = (
            gg.ggplot(
                gene_embedding_df.dropna(subset=["map_facet_label"]),
                gg.aes(x="umap_0", y="umap_1"),
            )
            + gg.geom_point(
                gg.aes(color="mAP", shape="p < 0.05"), size=2, stroke=0, alpha=0.5
            )
            + gg.facet_wrap("~map_facet_label")
            + gg.theme_bw()
            + gg.xlab("UMAP 0")
            + gg.ylab("UMAP 1")
            + gg.theme(strip_background=gg.element_rect(colour="black", fill="#fdfff4"))
            + gg.scale_color_gradient(low="blue", high="red", limits=[0, 1])
            + gg.scale_shape_manual(values={False: "s", True: "o"})
        )
        out_prefix = f"{gse_id}_{gene}_singlecell_umap_map"
        save_plot(gene_umap_gg, out_prefix, gene_umap_dir, dpi=500, height=5, width=6)


def main():
    set_plotting_style(font_size=5, linewidth=0.35)

    gse_id = "GSE132080"
    data_dir = Path("inputs")
    results_dir = Path("outputs")
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    colormap = make_colormap()
    scatterplot_size = (1.5, 1.25)
    faceted_size = (1.86, 1.5)

    bulk_df = load_bulk_map(gse_id, results_dir)
    create_bulk_scatter_plots(
        bulk_df,
        figures_dir,
        gse_id,
        colormap=colormap,
        scatterplot_size=scatterplot_size,
        faceted_size=faceted_size,
    )

    sc_df, _ = load_single_cell_data(gse_id, data_dir, results_dir, bulk_df)
    create_singlecell_global_plot(
        sc_df, figures_dir, gse_id, colormap=colormap, fig_size=scatterplot_size
    )
    create_singlecell_facet_plot(
        sc_df, figures_dir, gse_id, colormap=colormap, fig_size=faceted_size
    )
    plot_individual_umaps(sc_df, figures_dir, gse_id)


if __name__ == "__main__":
    main()
