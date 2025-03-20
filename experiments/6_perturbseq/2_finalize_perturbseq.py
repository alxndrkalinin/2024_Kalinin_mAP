# 2. Finalize the single cell perturbseq dataset
#
# The output of the Seurat pipeline in 1_process_perturbseq.ipynb is not easily compatible with our downstream tasks. This notebook finalizes the input perturbseq (CRISPRi) dataset.
#
# We use `pycytominer` to form bulk profiles from the single cell readouts.

from pathlib import Path

import pandas as pd
from pycytominer import aggregate


def load_gene_expression(gse_id: str, data_dir: Path) -> tuple[pd.DataFrame, list]:
    """
    Load and process the gene expression data from the processed matrix file.

    Returns:
        gene_exp_df: DataFrame with gene expression data.
        gene_features: Sorted list of measured gene features.
    """
    gene_exp_file = data_dir / f"{gse_id}_processed_matrix.tsv.gz"
    print(f"Loading gene expression data from {gene_exp_file}...")
    gene_exp_df = (
        pd.read_csv(gene_exp_file, sep="\t", index_col=0)
        .transpose()
        .reset_index()
        .rename({"index": "Metadata_barcode"}, axis="columns")
    )
    # Extract measured genes (all columns except Metadata_barcode)
    gene_features = gene_exp_df.columns.tolist()
    gene_features.remove("Metadata_barcode")

    # Create a new column "Metadata_sequence" using the barcode
    gene_exp_df = gene_exp_df.assign(
        Metadata_sequence=[x.split("-")[0] for x in gene_exp_df.Metadata_barcode]
    )
    gene_exp_df.columns.name = ""

    meta_features = ["Metadata_barcode", "Metadata_sequence"]
    gene_features = sorted(
        gene_exp_df.drop(meta_features, axis="columns").columns.tolist()
    )

    # Reorder columns so that meta features come first
    gene_exp_df = gene_exp_df.reindex(meta_features + gene_features, axis="columns")

    print("Gene expression data shape:", gene_exp_df.shape)
    print(gene_exp_df.head())
    return gene_exp_df, gene_features


def load_cell_identities(gse_id: str, data_dir: Path) -> pd.DataFrame:
    """
    Load and process cell identities data.

    Returns:
        DataFrame with cell identities.
    """
    identity_file = data_dir / f"{gse_id}_cell_identities.csv.gz"
    print(f"Loading cell identities from {identity_file}...")
    cell_id_df = pd.read_csv(identity_file, sep=",")

    # Prefix columns with "Metadata_"
    cell_id_df.columns = [f"Metadata_{x}" for x in cell_id_df.columns]
    # Create a column for gene identity from the guide identity column
    cell_id_df = cell_id_df.assign(
        Metadata_gene_identity=[
            str(x).split("_")[0] for x in cell_id_df.Metadata_guide_identity
        ]
    )

    print("Cell identities shape:", cell_id_df.shape)
    print(cell_id_df.head())
    return cell_id_df


def merge_single_cell_data(
    cell_id_df: pd.DataFrame, gene_exp_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge cell identity data with gene expression data.

    Returns:
        sc_df: DataFrame containing merged single cell data.
    """
    print("Merging cell identities with gene expression data...")
    sc_df = cell_id_df.merge(
        gene_exp_df,
        how="right",
        right_on="Metadata_barcode",
        left_on="Metadata_cell_barcode",
    )
    sc_df = sc_df.reset_index().rename(
        {"index": "Metadata_cell_identity"}, axis="columns"
    )
    # Label each cell with a unique identifier
    sc_df.Metadata_cell_identity = [
        f"sc_profile_{x}" for x in sc_df.Metadata_cell_identity
    ]

    print("Merged single cell data shape:", sc_df.shape)
    print(sc_df.head())
    return sc_df


def write_single_cell_data(sc_df: pd.DataFrame, output_file: Path) -> None:
    """
    Write the finalized single cell data to disk.
    """
    print(f"Writing single cell data to {output_file} ...")
    sc_df.to_csv(
        output_file,
        index=False,
        sep="\t",
        compression={"method": "gzip", "mtime": 1},
    )
    print("Single cell data written.")


def calculate_bulk_data(sc_df: pd.DataFrame, gene_features: list) -> pd.DataFrame:
    """
    Aggregate single cell data into bulk profiles using the median.

    Returns:
        bulk_df: DataFrame containing the aggregated bulk data.
    """
    print("Aggregating single cell data into bulk profiles...")
    bulk_df = aggregate(
        population_df=sc_df,
        strata=["Metadata_guide_identity"],
        features=gene_features,
        operation="median",
    )
    # Remove rows with null guide identity
    bulk_df = bulk_df[~bulk_df["Metadata_guide_identity"].isnull()]

    # Create a column for gene identity from guide identity and filter out unwanted entries
    bulk_df = bulk_df.assign(
        Metadata_gene_identity=[
            x.split("_")[0] for x in bulk_df.Metadata_guide_identity
        ]
    ).query("Metadata_gene_identity != '*'")

    # Reorder columns so that key identifiers come first
    bulk_df = bulk_df.reindex(
        ["Metadata_guide_identity", "Metadata_gene_identity"] + gene_features,
        axis="columns",
    )

    print("Bulk data shape:", bulk_df.shape)
    print(bulk_df.head())
    return bulk_df


def write_bulk_data(bulk_df: pd.DataFrame, output_bulk_file: Path) -> None:
    """
    Write the bulk aggregated data to disk.
    """
    print(f"Writing bulk data to {output_bulk_file} ...")
    bulk_df.to_csv(
        output_bulk_file,
        index=False,
        sep="\t",
        compression={"method": "gzip", "mtime": 1},
    )
    print("Bulk data written.")


def main():
    gse_id = "GSE132080"
    inputs_dir = Path("inputs")
    outputs_dir = Path("outputs")

    output_file = outputs_dir / f"{gse_id}_final_analytical.tsv.gz"
    output_bulk_file = outputs_dir / f"{gse_id}_bulk_final_analytical.tsv.gz"

    # Load gene expression and cell identities
    gene_exp_df, gene_features = load_gene_expression(gse_id, outputs_dir)
    cell_id_df = load_cell_identities(gse_id, inputs_dir)

    # Merge to create the single cell dataset
    sc_df = merge_single_cell_data(cell_id_df, gene_exp_df)

    # Write single cell data to disk
    write_single_cell_data(sc_df, output_file)

    # Aggregate single cell data to create bulk profiles
    bulk_df = calculate_bulk_data(sc_df, gene_features)
    # Optionally, you can inspect or query bulk_df here, e.g.,
    # print(bulk_df.query("Metadata_gene_identity == 'neg'"))

    # Write bulk data to disk
    write_bulk_data(bulk_df, output_bulk_file)


if __name__ == "__main__":
    main()
