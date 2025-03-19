#!/usr/bin/env python3
# coding: utf-8

# Preprocess Cell Health profiles, standardize or MAD robustize, and feature select.
#
# We use the MAD robustize in downstream analyses, but here we compare the impact
# of normalization strategy on plate and well position effect assesement (Figure 3).

from pathlib import Path

import pandas as pd
from tqdm import tqdm
from pycytominer import normalize, feature_select


def load_profile(plate, input_dir):
    file_path = input_dir / f"{plate}_augmented.csv.gz"
    return pd.read_csv(file_path)


def reformat_metadata(df):
    df = df.rename(
        {
            "Image_Metadata_Plate": "Metadata_Plate",
            "Image_Metadata_Well": "Metadata_Well",
        },
        axis="columns",
    ).drop(["Metadata_broad_sample"], axis="columns")
    metadata_cols = [col for col in df.columns if col.startswith("Metadata")]
    meta_df = df[metadata_cols]
    return pd.concat([meta_df, df.drop(columns=metadata_cols)], axis="columns")


def process_plate(plate, meta_features, method, input_dir, output_dir):
    df = load_profile(plate, input_dir)
    norm_file = output_dir / f"{plate}_wholeplate_normalized.csv.gz"
    normalize(
        profiles=df,
        features="infer",
        meta_features=meta_features,
        samples="all",
        method=method,
        output_file=norm_file,
        compression_options={"method": "gzip", "mtime": 1},
    )
    return norm_file


def merge_profiles(file_pattern, directory):
    files = sorted(directory.glob(file_pattern))
    merged = pd.concat([pd.read_csv(f) for f in files], sort=True)
    return reformat_metadata(merged)


def perform_feature_selection(df, ops, na_cutoff=0):
    df = feature_select(profiles=df, operation=ops, na_cutoff=na_cutoff)
    costes_cols = [col for col in df.columns if "costes" in col.lower()]
    return df.drop(columns=costes_cols)


def process_normalized_profiles(
    plates, meta_features, input_dir, output_dir, feature_select_ops
):
    for plate in tqdm(plates, desc="MAD robustizing profiles"):
        process_plate(
            plate,
            meta_features,
            method="mad_robustize",
            input_dir=input_dir,
            output_dir=output_dir,
        )
    merged_df = merge_profiles("*_normalized.csv.gz", output_dir)
    selected_df = perform_feature_selection(merged_df, feature_select_ops)
    out_file = (
        output_dir
        / "cell_health_profiles_merged_wholeplate_normalized_featureselected.tsv.gz"
    )
    selected_df.to_csv(out_file, index=False, sep="\t")


def process_standardized_profiles(
    plates, meta_features, input_dir, output_dir, feature_select_ops
):
    raw_df = pd.concat(
        [
            load_profile(plate, input_dir)
            for plate in tqdm(plates, desc="Loading raw profiles")
        ],
        sort=True,
    )
    raw_df = reformat_metadata(raw_df).reset_index(drop=True)
    meta_cols = [col for col in raw_df.columns if col.startswith("Metadata")]
    standardized_df = normalize(
        profiles=raw_df,
        features="infer",
        meta_features=meta_cols,
        samples="all",
        method="standardize",
    )
    selected_df = perform_feature_selection(standardized_df, feature_select_ops)
    out_file = (
        output_dir / "cell_health_profiles_merged_standardized_featureselected.tsv.gz"
    )
    selected_df.to_csv(out_file, index=False, sep="\t")


def main():
    input_dir = Path("inputs")
    output_dir = Path("outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    plates = [
        "SQ00014610",
        "SQ00014611",
        "SQ00014612",
        "SQ00014613",
        "SQ00014614",
        "SQ00014615",
        "SQ00014616",
        "SQ00014617",
        "SQ00014618",
    ]
    meta_features = [
        "Image_Metadata_Plate",
        "Image_Metadata_Well",
        "Metadata_WellRow",
        "Metadata_WellCol",
        "Metadata_gene_name",
        "Metadata_pert_name",
        "Metadata_broad_sample",
        "Metadata_cell_line",
    ]
    feature_select_ops = [
        "variance_threshold",
        "correlation_threshold",
        "drop_na_columns",
        "blocklist",
        "drop_outliers",
    ]

    process_normalized_profiles(
        plates, meta_features, input_dir, output_dir, feature_select_ops
    )
    process_standardized_profiles(
        plates, meta_features, input_dir, output_dir, feature_select_ops
    )


if __name__ == "__main__":
    main()
