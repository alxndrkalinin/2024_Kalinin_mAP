#!/usr/bin/env python
# coding: utf-8

# Assessing phenotypic consistency of CORUM complexes in the cpg0016[orf] dataset
#
# Based on the following notebook:
# https://github.com/jump-cellpainting/2024_Chandrasekaran_Morphmap/blob/96eff01ead6ea52dad11e9b13eed2f440ec1a8b4/03.retrieve-annotations/5.0-phenotypic-consistency-corum-complex.ipynb
#
# Chandrasekaran, S. N. et al. Morphological map of under- and over-expression of genes in human cells. bioRxiv 2024.12.02.624527 (2024)

from pathlib import Path

import pandas as pd
from copairs.map import mean_average_precision, multilabel
from orf_utils import consensus, get_metadata, get_featuredata


def load_orf_annotations(annotation_col: str, multi_label_col: str) -> pd.DataFrame:
    """Load and process the ORF annotations from the metadata file."""
    df = (
        pd.read_csv("inputs/orf_metadata.tsv.gz", sep="\t")[
            [
                "Metadata_JCP2022",
                annotation_col,
                "Metadata_pert_type",
                "Metadata_NCBI_Gene_ID",
            ]
        ]
        .dropna()
        .assign(col=lambda x: list(x[annotation_col].str.split("|")))
        .rename(columns={"col": multi_label_col})
    )
    print("ORF annotation dataframe shape:", df.shape)
    return df


def load_profiles(profiles: str) -> pd.DataFrame:
    """Load the ORF profiles from the parquet file."""
    parquet_file = f"inputs/profiles_{profiles}.parquet"
    df = pd.read_parquet(parquet_file)
    print("Profiles dataframe shape:", df.shape)
    return df


def load_phenotypic_activity(profiles: str, dir_path: Path) -> pd.DataFrame:
    """Load the phenotypic activity data and filter by below_corrected_p."""
    activity_file = dir_path / f"{profiles}_activity.parquet"
    df = pd.read_parquet(activity_file).query("below_corrected_p==True")
    print("Phenotypic activity dataframe shape:", df.shape)
    return df


def merge_data(
    profiles_df: pd.DataFrame,
    activity_df: pd.DataFrame,
    orf_annotation_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge profiles with phenotypic activity and ORF annotations."""
    df = profiles_df.merge(
        activity_df[["Metadata_JCP2022"]], on="Metadata_JCP2022", how="inner"
    )
    print("After merging phenotypic activity, shape:", df.shape)
    df = df.merge(orf_annotation_df, on="Metadata_JCP2022", how="inner")
    print("After merging ORF annotations, shape:", df.shape)
    return df


def compute_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Compute consensus by 'Metadata_NCBI_Gene_ID' and remove duplicate CORUM annotations."""
    consensus_df = consensus(df, "Metadata_NCBI_Gene_ID")
    # Remove duplicate entries in the CORUM annotations
    consensus_df.loc[:, "Metadata_corum_complex_list"] = (
        consensus_df.Metadata_corum_complex_list.apply(lambda x: list(set(x)))
    )
    print("Consensus dataframe shape:", consensus_df.shape)
    return consensus_df


def compute_map(
    consensus_df: pd.DataFrame,
    multi_label_col: str,
    batch_size: int,
    null_size: int,
    fdr: float,
) -> pd.DataFrame:
    """Compute multilabel average precision and mean average precision for CORUM complexes."""
    metadata_df = get_metadata(consensus_df)
    feature_df = get_featuredata(consensus_df)
    feature_values = feature_df.values

    pos_sameby = [multi_label_col]
    pos_diffby = []
    neg_sameby = []
    neg_diffby = [multi_label_col]

    result = multilabel.average_precision(
        metadata_df,
        feature_values,
        pos_sameby,
        pos_diffby,
        neg_sameby,
        neg_diffby,
        batch_size=batch_size,
        multilabel_col=multi_label_col,
    )
    print("Multilabel average precision result shape:", result.shape)

    agg_result = mean_average_precision(
        result, pos_sameby, null_size, threshold=fdr, seed=12527
    )
    print("Aggregated mAP result shape:", agg_result.shape)
    return agg_result


def main():
    profiles = "wellpos_cc_var_mad_outlier_featselect_sphering_harmony"
    retrieval_label = "corum_complex"
    multi_label_col = "Metadata_corum_complex_list"
    annotation_col = "Metadata_complexname"
    batch_size = 20000
    null_size = 20000
    fdr = 0.05
    output_dir = Path("outputs")

    orf_annotation_df = load_orf_annotations(annotation_col, multi_label_col)
    profiles_df = load_profiles(profiles)
    activity_df = load_phenotypic_activity(profiles, output_dir)
    merged_df = merge_data(profiles_df, activity_df, orf_annotation_df)
    consensus_df = compute_consensus(merged_df)

    # compute mAP for CORUM complexes
    agg_result = compute_map(consensus_df, multi_label_col, batch_size, null_size, fdr)

    output_file = output_dir / f"{retrieval_label}_consistency.parquet"
    agg_result.to_parquet(output_file, index=False)
    print(f"Phenotypic consistency mAP results saved to {output_file}")


if __name__ == "__main__":
    main()
