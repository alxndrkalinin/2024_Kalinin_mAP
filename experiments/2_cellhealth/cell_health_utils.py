import numpy as np
import plotnine as gg
from scipy.stats import combine_pvalues


def get_control_barcodes(profiles):
    # Define cell health constants
    barcode_col = "Metadata_pert_name"
    gene_col = "Metadata_gene_name"

    replicate_group = {"profile_col": barcode_col, "replicate_group_col": gene_col}

    control_group_cut = ["Chr2", "Luc", "LacZ"]
    control_group_pert = ["EMPTY"]

    control_barcodes_cut = (
        profiles.loc[
            profiles[replicate_group["replicate_group_col"]].isin(control_group_cut),
            replicate_group["profile_col"],
        ]
        .unique()
        .tolist()
    )

    control_barcodes_pert = (
        profiles.loc[
            profiles[replicate_group["replicate_group_col"]].isin(control_group_pert),
            replicate_group["profile_col"],
        ]
        .unique()
        .tolist()
    )

    control_barcodes = {
        "cutting_control": control_barcodes_cut,
        "perturbation_control": control_barcodes_pert,
    }

    return control_barcodes


def stouffer_method(p_values):
    _, combined_p_value = combine_pvalues(p_values, method="stouffer")
    return combined_p_value


def subset_6_replicates(profiles, cell_line=None, control="cutting_control"):
    control_barcodes = get_control_barcodes(profiles)
    if cell_line is not None:
        cell_line_df = profiles.query("Metadata_cell_line == @cell_line")
        n_replicates = 6
    else:
        cell_line_df = profiles
        n_replicates = 6 * profiles["Metadata_cell_line"].nunique()
    perts_6_wells = cell_line_df.groupby("Metadata_pert_name")["Metadata_Well"].count()
    perts_6_wells = perts_6_wells[perts_6_wells == n_replicates].index.tolist()
    df_6wells = cell_line_df[
        cell_line_df.Metadata_pert_name.isin(perts_6_wells)
        | cell_line_df.Metadata_pert_name.isin(control_barcodes[control])
    ].copy()
    df_6wells = df_6wells.reset_index(drop=True)
    df_6wells["Metadata_control_index"] = np.where(
        df_6wells.Metadata_pert_name.isin(control_barcodes[control]),
        df_6wells.index,
        -1,
    )
    df_6wells["Metadata_is_control"] = np.where(
        df_6wells.Metadata_pert_name.isin(control_barcodes[control]), 1, 0
    )
    return df_6wells


def get_cell_line_colors():
    return {"A549": "#1b9e77", "HCC44": "#d95f02", "ES2": "#7570b3"}


def plot_ap_per_label(
    ap_df, label_col, cell_line_colors, width=8, height=4, ncol=5, tick_length=1.0
):
    gg_ap = (
        gg.ggplot(ap_df, gg.aes(x="AP", y="-log10(AP p-value)", fill="Cell type"))
        + gg.geom_jitter(
            data=ap_df[~ap_df["p < 0.05"]], shape="s", size=1, color="white", stroke=0.1
        )
        + gg.geom_jitter(
            data=ap_df[ap_df["p < 0.05"]], shape="o", size=1, color="white", stroke=0.1
        )
        + gg.scale_fill_manual(values=cell_line_colors)
        + gg.geom_hline(
            yintercept=-np.log10(0.05), linetype="dashed", color="grey", size=0.2
        )
        + gg.theme_bw()
        + gg.xlab("AP")
        + gg.ylab("-log10(AP p-value)")
        + gg.facet_wrap(f"~{label_col}", ncol=ncol)
        + gg.scale_x_continuous(
            breaks=[0, 0.5, 1], labels=["0", "0.5", "1"], minor_breaks=[]
        )
        + gg.scale_y_continuous(
            breaks=[0, 2, 4], labels=["0", "2", "4"], minor_breaks=[]
        )
        + gg.theme(
            strip_background=gg.element_rect(
                colour="black", fill="#fdfff4", linewidth=0.2
            ),
            strip_text=gg.element_text(
                margin={"t": 0, "b": 0}, size=5, family="Open Sans"
            ),
            panel_border=gg.element_rect(colour="grey", fill=None, linewidth=0.2),
            panel_background=gg.element_rect(fill="white", colour=None),
            plot_background=gg.element_rect(fill="white", colour=None),
            panel_grid_major=gg.element_blank(),
            panel_grid_minor=gg.element_blank(),
            axis_text=gg.element_text(family="Open Sans", size=5),
            text=gg.element_text(family="Open Sans", size=5),
            axis_title=gg.element_text(family="Open Sans", size=5),
            figure_size=(width, height),
            legend_position="none",
            axis_ticks=gg.element_line(linewidth=0.2),
            axis_line=gg.element_line(linewidth=0.2),
        )
    )

    fig = gg_ap.draw()
    fig.set_size_inches(width, height)
    for ax in fig.axes:
        ax.tick_params(length=tick_length)

    return fig
