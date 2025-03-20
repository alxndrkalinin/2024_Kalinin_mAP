#!/usr/bin/env python3
# coding: utf-8

# Assessment of gene phenotypic consistency (Figs 3E & 3F)
#
# Calculate mAP to assess phenotypic consistency of each gene (retrieval its guides against other genes' guides).

from pathlib import Path

import numpy as np
import pandas as pd
import plotnine as gg
from pycytominer import aggregate

from map_utils.map import calculate_map
from map_utils.plot import plot_map_scatter, set_plotting_style
from cell_health_utils import subset_6_replicates, get_cell_line_colors


def save_plot(plot_obj, base_name, output_dir, dpi=500, width=8, height=4):
    output_dir.mkdir(parents=True, exist_ok=True)
    svg_path = output_dir / f"{base_name}.svg"
    png_path = output_dir / f"{base_name}.png"
    if hasattr(plot_obj, "save"):
        plot_obj.save(
            str(svg_path), dpi=dpi, height=height, width=width, bbox_inches="tight"
        )
        plot_obj.save(
            str(png_path), dpi=dpi, height=height, width=width, bbox_inches="tight"
        )
    elif hasattr(plot_obj, "savefig"):
        plot_obj.savefig(str(svg_path), dpi=dpi, bbox_inches="tight")
        plot_obj.savefig(str(png_path), dpi=dpi, bbox_inches="tight")
    else:
        raise TypeError("Plot object does not have a save or savefig method.")


def load_data(profiles_path, phen_activity_path):
    profiles = pd.read_csv(profiles_path, sep="\t")
    phen_activity_results = pd.read_csv(phen_activity_path)
    return profiles, phen_activity_results


def get_replicable_perts(phen_activity_results):
    replicable_perts = (
        phen_activity_results[phen_activity_results["p < 0.05"]]
        .groupby("Cell type")["Metadata_pert_name"]
        .unique()
    )
    return replicable_perts


def get_aggregated_data(profiles, cell_lines):
    df_list = [subset_6_replicates(profiles, cell) for cell in cell_lines]
    ch_df_6well = pd.concat(df_list, axis=0).drop_duplicates()
    ch_df_agg = aggregate(
        ch_df_6well,
        strata=["Metadata_pert_name", "Metadata_gene_name", "Metadata_cell_line"],
        features="infer",
    )
    assert ch_df_agg.index.is_unique
    return ch_df_agg


def get_replicable_genes(ch_df_agg, replicable_perts):
    replicable_genes = {}
    for cell in ch_df_agg.Metadata_cell_line.unique():
        vc = ch_df_agg[
            (ch_df_agg["Metadata_cell_line"] == cell)
            & (ch_df_agg["Metadata_pert_name"].isin(replicable_perts[cell]))
        ].Metadata_gene_name.value_counts()
        replicable_genes[cell] = vc[vc > 1].index.tolist()
    return replicable_genes


def compute_map_metrics(ch_df_agg, replicable_genes, cell_lines):
    pair_config = {
        "pos_sameby": ["Metadata_gene_name"],
        "pos_diffby": [],
        "neg_sameby": [],
        "neg_diffby": ["Metadata_gene_name"],
    }
    map_config = {"null_size": 1000000, "groupby_columns": ["Metadata_gene_name"]}
    ap_results_list = []
    map_results_list = []
    for cell in ch_df_agg.Metadata_cell_line.unique():
        print(f"Processing cell type: {cell}")
        df = ch_df_agg.query(
            "Metadata_cell_line == @cell and Metadata_gene_name in @replicable_genes[@cell]"
        ).reset_index(drop=True)
        map_scores, ap_scores = calculate_map(
            df, pair_config, map_config, return_ap_scores=True
        )
        ap_scores["Cell type"] = cell
        map_scores["Cell type"] = cell
        ap_results_list.append(ap_scores)
        map_results_list.append(map_scores)
    ap_results = pd.concat(ap_results_list)
    ap_results["markers"] = np.where(ap_results["p_value"] < 0.05, "p<0.05", "p>=0.05")
    ap_results["const"] = 1
    ap_results.rename(columns={"average_precision": "AP"}, inplace=True)
    map_results = pd.concat(map_results_list)
    map_results["markers"] = np.where(
        map_results["p_value"] < 0.05, "p<0.05", "p>=0.05"
    )
    map_results["const"] = 1
    return ap_results, map_results


def apply_jitter(map_results, pvalue_jitter_strength=0.2, map_jitter_strength=0.05):
    results_jitter = map_results.copy()
    results_jitter["-log10(mAP p-value)"] = (
        results_jitter["-log10(mAP p-value)"]
        + np.random.uniform(
            -pvalue_jitter_strength, pvalue_jitter_strength, results_jitter.shape[0]
        )
    ).clip(0, 5)
    results_jitter["mAP"] = (
        results_jitter["mAP"]
        + np.random.uniform(
            -map_jitter_strength, map_jitter_strength, results_jitter.shape[0]
        )
    ).clip(0, 1)
    return results_jitter


def plot_fig3E(results_jitter, cell_line_colors, fig_dir):
    plot_obj = plot_map_scatter(
        results_jitter,
        "const",
        "",
        hue_col="Cell type",
        palette=cell_line_colors,
        row=None,
        move_legend="center left",
        pr_x=0.1,
        pr_y=0.75,
        l_x=1.15,
        l_y=0.8,
        figure="Fig3E",
        save_path=str(fig_dir),
    )
    save_plot(plot_obj, "Fig3E", fig_dir, dpi=500, width=8, height=4)


def plot_fig3F(ap_results, map_results, cell_line_colors, fig_dir):
    gene_counts = map_results[map_results["p < 0.05"]].Metadata_gene_name.value_counts()
    genes_consistent = gene_counts[gene_counts > 1].index.tolist()  # noqa
    consistent_results = ap_results.query("Metadata_gene_name in @genes_consistent")
    gene_gg = (
        gg.ggplot(
            consistent_results,
            gg.aes(x="AP", y="-log10(AP p-value)", color="Cell type"),
        )
        + gg.geom_jitter(
            data=consistent_results[~consistent_results["p < 0.05"]],
            shape="s",
            size=1.2,
        )
        + gg.geom_jitter(
            data=consistent_results[consistent_results["p < 0.05"]], shape="o", size=1.2
        )
        + gg.geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="grey")
        + gg.theme_bw()
        + gg.xlab("AP")
        + gg.ylab("-log10(AP p-value)")
        + gg.facet_wrap("~Metadata_gene_name", ncol=5)
        + gg.scale_color_manual(values=cell_line_colors)
        + gg.scale_x_continuous(breaks=[0, 0.5, 1.0], labels=["0", "0.5", "1"])
        + gg.theme(
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
            axis_text=gg.element_text(size=10),
            text=gg.element_text(family="Open Sans", size=14),
            axis_title=gg.element_text(family="Open Sans", size=14),
            legend_title=gg.element_text(margin={"b": 20}),
            strip_text=gg.element_text(size=10, family="Open Sans"),
            figure_size=(8, 4),
        )
    )
    save_plot(gene_gg, "Fig3F", fig_dir, dpi=500, width=8, height=4)


def main():
    set_plotting_style()
    cell_lines = ["A549", "ES2", "HCC44"]
    cell_line_colors = get_cell_line_colors()

    output_dir = Path("outputs")
    profiles_path = (
        output_dir
        / "cell_health_profiles_merged_wholeplate_normalized_featureselected.tsv.gz"
    )
    phen_activity_path = output_dir / "phenotypic_activity_map.csv"
    profiles, phen_activity_results = load_data(profiles_path, phen_activity_path)
    print(profiles.shape)

    replicable_perts = get_replicable_perts(phen_activity_results)
    ch_df_agg = get_aggregated_data(profiles, cell_lines)
    replicable_genes = get_replicable_genes(ch_df_agg, replicable_perts)

    ap_results, map_results = compute_map_metrics(
        ch_df_agg, replicable_genes, cell_lines
    )
    results_jitter = apply_jitter(map_results)

    fig_dir = Path("outputs/figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_fig3E(results_jitter, cell_line_colors, fig_dir)
    plot_fig3F(ap_results, map_results, cell_line_colors, fig_dir)


if __name__ == "__main__":
    main()
