#!/usr/bin/env python3
# coding: utf-8

# Phenotypic activity (Figs 3B, 3C)
#
# Calculate mAP to assess phenotypic activity of each perturbation,
#  i.e. retrieval its replicates against negative controls.

from pathlib import Path

import mygene
import numpy as np
import pandas as pd
import plotnine as gg

from map_utils.map import calculate_map
from map_utils.plot import plot_map_scatter, set_plotting_style
from cell_health_utils import subset_6_replicates, get_cell_line_colors


def save_plot(plot_obj, base_name, output_dir, dpi=500, width=8, height=4):
    svg_path = output_dir / f"{base_name}.svg"
    png_path = output_dir / f"{base_name}.png"
    # If the plot object has a 'save' method (e.g., plotnine), use it.
    if hasattr(plot_obj, "save"):
        plot_obj.save(
            str(svg_path), dpi=dpi, height=height, width=width, bbox_inches="tight"
        )
        plot_obj.save(
            str(png_path), dpi=dpi, height=height, width=width, bbox_inches="tight"
        )
    # Otherwise, assume it's a matplotlib figure.
    elif hasattr(plot_obj, "savefig"):
        plot_obj.savefig(str(svg_path), dpi=dpi, bbox_inches="tight")
        plot_obj.savefig(str(png_path), dpi=dpi, bbox_inches="tight")
    else:
        raise TypeError(
            "The provided plot object does not have a save or savefig method."
        )


def process_phenotypic_activity(profiles, cell_lines, results_dir):
    pair_config = {
        "pos_sameby": ["Metadata_pert_name", "Metadata_control_index"],
        "pos_diffby": [],
        "neg_sameby": [],
        "neg_diffby": [
            "Metadata_pert_name",
            "Metadata_control_index",
            "Metadata_is_control",
        ],
    }
    map_config = {"null_size": 1000000, "groupby_columns": ["Metadata_pert_name"]}
    ap_list, map_list = [], []
    for cell in cell_lines:
        df6 = subset_6_replicates(profiles, cell)
        print(f"Processing cell type: {cell}")
        map_scores, ap_scores = calculate_map(
            df6, pair_config, map_config, return_ap_scores=True
        )
        ap_scores["Cell type"] = cell
        map_scores["Cell type"] = cell
        ap_list.append(ap_scores)
        map_list.append(map_scores)
    ap_res = pd.concat(ap_list).reset_index(drop=True)
    ap_res["markers"] = np.where(ap_res["p_value"] < 0.05, "p<0.05", "p>=0.05")
    map_res = pd.concat(map_list).reset_index(drop=True)
    map_res["markers"] = np.where(map_res["p_value"] < 0.05, "p<0.05", "p>=0.05")
    ap_res.to_csv(results_dir / "phenotypic_activity_ap.csv", index=False)
    map_res.to_csv(results_dir / "phenotypic_activity_map.csv", index=False)
    return ap_res, map_res


def plot_map_figure(map_res, cell_line_colors, fig_dir):
    map_res["const"] = 1
    fig_map = plot_map_scatter(
        map_res,
        "const",
        "",
        hue_col="Cell type",
        palette=cell_line_colors,
        row=None,
        move_legend="center left",
        pr_x=0.53,
        pr_y=0.3,
        l_x=1.15,
        l_y=0.8,
        figure="Fig3B",
        save_path=str(fig_dir),
    )
    return fig_map


def plot_individual_guides(map_res, cell_line_colors, fig_dir):
    guide_counts = map_res[map_res["p < 0.05"]].Metadata_pert_name.value_counts()
    guides = sorted(guide_counts[guide_counts > 1].index)[:18]
    guide_res = map_res.query("Metadata_pert_name in @guides")
    gene_gg = (
        gg.ggplot(
            guide_res, gg.aes(x="mAP", y="-log10(mAP p-value)", color="Cell type")
        )
        + gg.geom_point(data=guide_res[~guide_res["p < 0.05"]], shape="s", size=1.5)
        + gg.geom_point(data=guide_res[guide_res["p < 0.05"]], shape="o", size=1.5)
        + gg.geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="grey")
        + gg.theme_bw()
        + gg.xlab("mAP")
        + gg.ylab("-log10(mAP p-value)")
        + gg.facet_wrap("~Metadata_pert_name", ncol=6)
        + gg.scale_color_manual(values=cell_line_colors)
        + gg.scale_x_continuous(breaks=[0, 0.5, 1], labels=["0", "0.5", "1"])
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
    save_plot(gene_gg, "Fig3C", fig_dir, dpi=500, width=8, height=4)
    return guides


def plot_guide_ap(ap_res, cell_line_colors, guides, fig_dir):
    ap_res.rename(columns={"average_precision": "AP"}, inplace=True)
    guide_ap = ap_res.query("Metadata_pert_name in @guides")
    gene_gg_ap = (
        gg.ggplot(guide_ap, gg.aes(x="AP", y="-log10(AP p-value)", color="Cell type"))
        + gg.geom_jitter(data=guide_ap[~guide_ap["p < 0.05"]], shape="s", size=1.2)
        + gg.geom_jitter(data=guide_ap[guide_ap["p < 0.05"]], shape="o", size=1.2)
        + gg.geom_hline(yintercept=-np.log10(0.05), linetype="dashed", color="grey")
        + gg.theme_bw()
        + gg.xlab("AP")
        + gg.ylab("-log10(AP p-value)")
        + gg.facet_wrap("~Metadata_pert_name", ncol=6)
        + gg.scale_color_manual(values=cell_line_colors)
        + gg.scale_x_continuous(breaks=[0, 0.5, 1], labels=["0", "0.5", "1"])
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
    save_plot(gene_gg_ap, "Fig3C_AP", fig_dir, dpi=500, width=8, height=4)


def compare_controls(profiles, map_res, cell_lines, fig_dir, results_dir):
    empty_list = []
    for cell in cell_lines:
        df_empty = subset_6_replicates(profiles, cell, control="perturbation_control")
        df_empty = df_empty.query("Metadata_cell_line == @cell").reset_index(drop=True)
        map_scores = calculate_map(
            df_empty,
            {
                "pos_sameby": ["Metadata_pert_name", "Metadata_control_index"],
                "pos_diffby": [],
                "neg_sameby": [],
                "neg_diffby": [
                    "Metadata_pert_name",
                    "Metadata_control_index",
                    "Metadata_is_control",
                ],
            },
            {"null_size": 1000000, "groupby_columns": ["Metadata_pert_name"]},
        )
        map_scores["cell_type"] = cell
        empty_list.append(map_scores)
    empty_res = pd.concat(empty_list).reset_index(drop=True)
    empty_res["markers"] = np.where(empty_res["p_value"] < 0.05, "p<0.05", "p>=0.05")
    empty_res.rename(columns={"cell_type": "Cell type"}, inplace=True)

    map_copy = map_res.copy()
    map_copy["control_barcodes"] = "cutting_control"
    empty_res["control_barcodes"] = "perturbation_control"
    guide_res = pd.concat([map_copy, empty_res]).reset_index(drop=True)
    guide_res = guide_res.merge(
        profiles[["Metadata_gene_name", "Metadata_pert_name"]].drop_duplicates(),
        on="Metadata_pert_name",
        how="left",
    )
    guide_res.to_csv(results_dir / "phenotypic_activity_empty.csv", index=False)

    control_df = guide_res.pivot_table(
        index=["Metadata_pert_name", "Metadata_gene_name", "Cell type"],
        columns="control_barcodes",
        values="-log10(mAP p-value)",
    ).reset_index()
    cutting_mean = guide_res.query("control_barcodes=='cutting_control'")["mAP"].mean()
    pert_mean = guide_res.query("control_barcodes=='perturbation_control'")[
        "mAP"
    ].mean()
    above_diag = (control_df["cutting_control"] > -np.log10(0.05)).mean()
    below_diag = (control_df["perturbation_control"] > -np.log10(0.05)).mean()

    comp_gg = (
        gg.ggplot(control_df, gg.aes(x="perturbation_control", y="cutting_control"))
        + gg.geom_point(gg.aes(fill="Cell type"), size=3, stroke=0.0, alpha=0.5)
        + gg.scale_fill_manual(name="Cell type", values=get_cell_line_colors())
        + gg.theme_bw()
        + gg.xlab("-log10(mAP p-value), perturbation control (empty)")
        + gg.ylab("-log10(mAP p-value), cutting control")
        + gg.geom_abline(intercept=0, slope=1, linetype="dashed", color="red")
        + gg.coord_fixed()
        + gg.annotate(
            "text",
            x=control_df["perturbation_control"].min(),
            y=control_df["cutting_control"].max(),
            label=f"Mean mAP: {cutting_mean:.2f}\nRetrieved: {below_diag:.0%}",
            ha="left",
            va="top",
        )
        + gg.annotate(
            "text",
            x=control_df["perturbation_control"].max(),
            y=control_df["cutting_control"].min(),
            label=f"Mean mAP: {pert_mean:.2f}\nRetrieved: {above_diag:.0%}",
            ha="right",
            va="bottom",
        )
        + gg.theme(legend_key=gg.element_rect(color="white"))
    )
    save_plot(comp_gg, "control_comparison", fig_dir, dpi=500, width=4, height=3.5)


def plot_phenotypic_vs_ceres(profiles, guide_res_file, cell_lines, fig_dir):
    mg = mygene.MyGeneInfo()
    guide_res = pd.read_csv(guide_res_file)
    result = mg.querymany(
        guide_res["Metadata_gene_name"].unique().tolist(),
        scopes="symbol",
        species="human",
        fields="entrezgene,symbol,ensembl.gene,",
        as_dataframe=True,
    )
    # Reset index to get a "query" column for merging.
    ncbi_df = result.reset_index().drop_duplicates(subset="_id")[["query", "_id"]]
    guide_res = guide_res.merge(
        ncbi_df, left_on="Metadata_gene_name", right_on="query", how="left"
    )

    ceres_dir = Path("inputs")
    ceres_df = pd.read_csv(ceres_dir / "ceres.csv", index_col=0)
    depmap_df = pd.read_csv(ceres_dir / "depmap_sample_info.csv", index_col=0)
    ncbi_ids = [col.split(" ")[1].strip("()") for col in ceres_df.columns]
    ceres_df.columns = ncbi_ids
    ceres_df = depmap_df.merge(ceres_df, left_index=True, right_index=True, how="right")
    assert all(
        cell in ceres_df["stripped_cell_line_name"].tolist() for cell in cell_lines
    )

    cols = ["stripped_cell_line_name"] + guide_res["_id"].dropna().unique().tolist()
    ceres_sub = (
        ceres_df.query("stripped_cell_line_name in @cell_lines")
        .loc[:, cols]
        .reset_index()
        .melt(
            id_vars=["DepMap_ID", "stripped_cell_line_name"],
            var_name="_id",
            value_name="ceres_score",
        )
    )
    merged_df = guide_res.merge(
        ceres_sub,
        left_on=["_id", "Cell type"],
        right_on=["_id", "stripped_cell_line_name"],
        how="left",
    )
    map_ceres_gg = (
        gg.ggplot(
            merged_df, gg.aes(x="-log10(mAP p-value)", y="ceres_score", size="mAP")
        )
        + gg.geom_point(gg.aes(fill="Cell type"), stroke=0.0, alpha=0.3)
        + gg.scale_fill_manual(name="Cell type", values=get_cell_line_colors())
        + gg.theme_bw()
        + gg.xlab("-log10(mAP p-value)")
        + gg.ylab("CERES")
        + gg.facet_wrap("~control_barcodes", ncol=2)
        + gg.theme(
            strip_background=gg.element_rect(color="black", fill="#fdfff4"),
            legend_key=gg.element_rect(color="white"),
        )
    )
    save_plot(
        map_ceres_gg, "map_ceres_comparison", fig_dir, dpi=500, width=6, height=3.5
    )


def main():
    set_plotting_style()
    cell_line_colors = get_cell_line_colors()

    profiles_file = Path(
        "outputs/cell_health_profiles_merged_wholeplate_normalized_featureselected.tsv.gz"
    )
    profiles = pd.read_csv(profiles_file, sep="\t")
    print(profiles.shape)

    results_dir = Path("outputs")
    results_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = results_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    cell_lines = ["A549", "ES2", "HCC44"]
    ap_res, map_res = process_phenotypic_activity(profiles, cell_lines, results_dir)
    fig_map = plot_map_figure(map_res, cell_line_colors, fig_dir)
    save_plot(fig_map, "Fig3B", fig_dir, dpi=500, width=8, height=4)

    guides = plot_individual_guides(map_res, cell_line_colors, fig_dir)
    plot_guide_ap(ap_res, cell_line_colors, guides, fig_dir)

    compare_controls(profiles, map_res, cell_lines, fig_dir, results_dir)

    guide_csv = results_dir / "phenotypic_activity_empty.csv"
    plot_phenotypic_vs_ceres(profiles, guide_csv, cell_lines, fig_dir)


if __name__ == "__main__":
    main()
