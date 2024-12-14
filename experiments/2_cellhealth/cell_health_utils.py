import numpy as np
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


def get_6well_plate_pairs(profiles, cell_line, control="cutting_control"):
    control_barcodes = get_control_barcodes(profiles)
    cell_line_df = profiles.query("Metadata_cell_line == @cell_line")
    perts_6_wells = cell_line_df.groupby("Metadata_pert_name")["Metadata_Well"].count()
    perts_6_wells = perts_6_wells[perts_6_wells == 6].index.tolist()
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
    return df_6wells


def get_cell_line_colors():
    return {"A549": "#1b9e77", "HCC44": "#d95f02", "ES2": "#7570b3"}
