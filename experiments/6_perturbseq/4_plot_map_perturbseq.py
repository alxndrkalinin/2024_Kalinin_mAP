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

from map_utils.plot import set_plotting_style


def save_plot(
    plot_obj, filename_prefix: str, output_dir: Path, dpi=500, height=5, width=6
):
    """
    Save a plotnine plot to both PNG and SVG formats.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    png_file = output_dir / f"{filename_prefix}.png"
    svg_file = output_dir / f"{filename_prefix}.svg"
    plot_obj.save(
        filename=str(png_file), dpi=dpi, height=height, width=width, bbox_inches="tight"
    )
    plot_obj.save(
        filename=str(svg_file), dpi=dpi, height=height, width=width, bbox_inches="tight"
    )
    print(f"Saved plot as:\n  {png_file}\n  {svg_file}")


def load_bulk_map(gse_id: str, results_dir: Path) -> pd.DataFrame:
    """
    Load bulk mAP results from file and process the DataFrame.
    """
    bulk_file = results_dir / f"{gse_id}_map.tsv"
    bulk_df = pd.read_csv(bulk_file, sep="\t")
    # Set gene as a categorical variable
    bulk_df["gene"] = pd.Categorical(bulk_df.gene, categories=bulk_df.gene.unique())
    print("Bulk map data shape:", bulk_df.shape)
    print(bulk_df.head(2))
    return bulk_df


def create_bulk_scatter_plots(bulk_df: pd.DataFrame, output_dir: Path, gse_id: str):
    """
    Create and save global and facet scatter plots for bulk data.
    """
    # Define color map
    palette = sns.color_palette().as_hex()
    color_map = {False: palette[0], True: palette[1]}

    # Global scatter plot (with linear regression smoothing)
    global_gg = (
        gg.ggplot(
            bulk_df, gg.aes(x="relative_activity_day5", y="mAP", color="p < 0.05")
        )
        + gg.geom_point(alpha=1.0, size=1)
        + gg.geom_smooth(method="lm", color="black")
        + gg.theme_bw()
        + gg.xlab("Mismatched guide RNA activity relative to perfect match")
        + gg.ylab("mAP")
        + gg.scale_color_manual(values=color_map)
    )
    save_plot(
        global_gg,
        f"{gse_id}_crispri_map_relative_activity_global_bulk",
        output_dir,
        dpi=500,
        height=3.5,
        width=4,
    )

    # Faceted scatter plot by gene
    gene_gg = (
        gg.ggplot(
            bulk_df, gg.aes(x="relative_activity_day5", y="mAP", color="p < 0.05")
        )
        + gg.geom_point(size=1.0)
        + gg.theme_bw()
        + gg.xlab("Mismatched guide RNA activity relative to perfect match")
        + gg.ylab("mAP")
        + gg.scale_color_manual(values=color_map)
        + gg.guides(color=gg.guide_colorbar(title="mAP"))
        + gg.facet_wrap("~gene")
        + gg.theme(
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            axis_text=gg.element_text(size=10),
        )
    )
    save_plot(
        gene_gg,
        f"{gse_id}_crispri_map_relative_activity_facet_bulk",
        output_dir,
        dpi=500,
        height=5,
        width=6,
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


def create_singlecell_global_plot(sc_df: pd.DataFrame, output_dir: Path, gse_id: str):
    """
    Create and save a global scatter plot (density and points) for single cell data.
    """
    # Define color palette (reuse from bulk)
    palette = sns.color_palette().as_hex()
    color_map = {False: palette[0], True: palette[1]}

    global_gg = (
        gg.ggplot(
            sc_df.dropna(subset=["Metadata_gene_identity"]),
            gg.aes(x="relative_activity_day5", y="mAP", color="p < 0.05"),
        )
        + gg.geom_density_2d()
        + gg.geom_point(size=0.2, alpha=0.1)
        + gg.theme_bw()
        + gg.xlab("Mismatched guide RNA activity relative to perfect match")
        + gg.ylab("mAP")
        + gg.scale_color_manual(values=color_map)
        + gg.theme(
            text=gg.element_text(family="Open Sans", size=14),
            axis_title=gg.element_text(family="Open Sans", size=14),
            legend_title=gg.element_text(margin={"b": 20}),
        )
    )
    save_plot(
        global_gg, f"{gse_id}_singlecell_global", output_dir, dpi=500, height=5, width=6
    )


def create_singlecell_facet_plot(sc_df: pd.DataFrame, output_dir: Path, gse_id: str):
    """
    Create and save a faceted scatter plot for single cell data by gene.
    """
    palette = sns.color_palette().as_hex()
    color_map = {False: palette[0], True: palette[1]}

    facet_gg = (
        gg.ggplot(
            sc_df.dropna(subset=["gene"]),
            gg.aes(x="relative_activity_day5", y="mAP", color="p < 0.05"),
        )
        + gg.geom_density_2d()
        + gg.geom_point(alpha=0.05, size=0.1)
        + gg.theme_bw()
        + gg.xlab("Mismatched guide RNA activity relative to perfect match")
        + gg.ylab("mAP")
        + gg.scale_color_manual(values=color_map)
        + gg.facet_wrap("~gene")
        + gg.theme(
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            axis_text=gg.element_text(size=10),
            text=gg.element_text(family="Open Sans", size=14),
            axis_title=gg.element_text(family="Open Sans", size=14),
            legend_title=gg.element_text(margin={"b": 20}),
            strip_text=gg.element_text(size=10, family="Open Sans"),
        )
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
    set_plotting_style()

    gse_id = "GSE132080"
    data_dir = Path("inputs")
    results_dir = Path("outputs")
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    bulk_df = load_bulk_map(gse_id, results_dir)
    create_bulk_scatter_plots(bulk_df, figures_dir, gse_id)

    sc_df, _ = load_single_cell_data(gse_id, data_dir, results_dir, bulk_df)
    create_singlecell_global_plot(sc_df, figures_dir, gse_id)
    create_singlecell_facet_plot(sc_df, figures_dir, gse_id)
    plot_individual_umaps(sc_df, figures_dir, gse_id)


if __name__ == "__main__":
    main()
