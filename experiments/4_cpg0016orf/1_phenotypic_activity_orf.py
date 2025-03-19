#!/usr/bin/env python
# coding: utf-8

# Assessing phenotypic activity of ORF perturbations in the cpg0016[orf] dataset
#
# Based on the following notebook:
# https://github.com/jump-cellpainting/2024_Chandrasekaran_Morphmap/blob/96eff01ead6ea52dad11e9b13eed2f440ec1a8b4/03.retrieve-annotations/0.0.phenotypic-activity-orf.ipynb
#
# Chandrasekaran, S. N. et al. Morphological map of under- and over-expression of genes in human cells. bioRxiv 2024.12.02.624527 (2024)

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from copairs.map import average_precision, mean_average_precision
from orf_utils import (
    get_featuredata,
    get_metadata,
    remove_empty_wells,
    remove_low_infection_efficiency_wells,
    remove_nan_features,
)

warnings.simplefilter(action="ignore", category=FutureWarning)


def load_orf_data(operations: str) -> pd.DataFrame:
    """Load the ORF parquet file."""
    orf_df = pd.read_parquet(f"inputs/profiles_{operations}.parquet")
    print("Loaded ORF dataframe shape:", orf_df.shape)
    return orf_df


def add_annotations(orf_df: pd.DataFrame) -> pd.DataFrame:
    """Add annotations from metadata files to the ORF dataframe."""
    orf_metdata_df = pd.read_csv("inputs/orf_metadata.tsv.gz", sep="\t")
    compound_metadata_df = pd.read_csv(
        "inputs/compound.csv.gz", usecols=["Metadata_JCP2022"]
    ).assign(
        Metadata_pert_type=lambda x: np.where(
            x["Metadata_JCP2022"] == "JCP2022_999999", "empty", "poscon"
        )
    )
    metadata_df = pd.concat(
        [orf_metdata_df, compound_metadata_df], join="outer", ignore_index=True
    )
    orf_df = orf_df.merge(metadata_df, on="Metadata_JCP2022", how="inner")
    print("After merging annotations, ORF dataframe shape:", orf_df.shape)
    return orf_df


def filter_data(orf_df: pd.DataFrame) -> pd.DataFrame:
    """Apply various filters to clean the ORF dataframe."""
    # Remove empty wells
    orf_df = remove_empty_wells(orf_df)
    print("After removing empty wells, shape:", orf_df.shape)

    # Remove 'poscon' wells
    orf_df = orf_df.query('Metadata_pert_type!="poscon"').reset_index(drop=True)
    print("After removing poscon wells, shape:", orf_df.shape)

    # Remove 'BAD CONSTRUCT' profiles
    orf_df = orf_df.query('Metadata_broad_sample!="BAD CONSTRUCT"').reset_index(
        drop=True
    )
    print("After removing BAD CONSTRUCT profiles, shape:", orf_df.shape)

    # Remove features with NaN values
    orf_df = remove_nan_features(orf_df)

    # Merge with platemap metadata and remove low infection efficiency wells
    platemap_df = pd.read_csv(
        "inputs/experiment-metadata.tsv",
        sep="\t",
        usecols=["Plate_Map_Name", "Assay_Plate_Barcode"],
    ).rename(
        columns={
            "Plate_Map_Name": "Metadata_plate_map_name",
            "Assay_Plate_Barcode": "Metadata_Plate",
        }
    )
    orf_df = orf_df.merge(platemap_df, on="Metadata_Plate", how="left")
    orf_df = remove_low_infection_efficiency_wells(orf_df)
    print("After removing low infection efficiency wells, shape:", orf_df.shape)
    return orf_df


def calculate_map(
    orf_df: pd.DataFrame, batch_size: int, null_size: int, fdr: float
) -> pd.DataFrame:
    """Calculate mean average precision (mAP) for each ORF perturbation."""
    # Add a column for negative control
    orf_df["Metadata_negcon"] = np.where(orf_df["Metadata_pert_type"] == "negcon", 1, 0)

    # Define grouping variables
    pos_sameby = ["Metadata_JCP2022"]
    pos_diffby = []
    neg_sameby = ["Metadata_Plate"]
    neg_diffby = ["Metadata_negcon"]

    # Retrieve metadata and feature data
    metadata_df = get_metadata(orf_df)
    feature_df = get_featuredata(orf_df)
    feature_values = feature_df.values

    # Calculate average precision
    result = average_precision(
        metadata_df,
        feature_values,
        pos_sameby,
        pos_diffby,
        neg_sameby,
        neg_diffby,
        batch_size=batch_size,
    )
    print("Average precision result shape:", result.shape)

    # Remove negative controls from result
    result = result.query('Metadata_pert_type!="negcon"').reset_index(drop=True)
    print("After removing negcon, result shape:", result.shape)

    # Calculate mean average precision
    agg_result = mean_average_precision(
        result, pos_sameby, null_size=null_size, threshold=fdr, seed=12527
    ).rename(columns={"average_precision": "mean_average_precision"})
    print("Aggregated mAP result shape:", agg_result.shape)
    return agg_result


def main():
    operations = "wellpos_cc_var_mad_outlier_featselect_sphering_harmony"
    batch_size = 20000
    null_size = 20000
    fdr = 0.1
    output_dir = Path("outputs")

    orf_df = load_orf_data(operations)
    orf_df = add_annotations(orf_df)
    orf_df = filter_data(orf_df)

    agg_result = calculate_map(orf_df, batch_size, null_size, fdr)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{operations}_activity.parquet"
    agg_result.to_parquet(output_file, index=False)
    print(f"Phenotypic activity mAP results saved to {output_file}")


if __name__ == "__main__":
    main()
