#!/usr/bin/env python3
# coding: utf-8

from pathlib import Path

import pandas as pd
from pycytominer import feature_select


def split_data(df: pd.DataFrame, dataset: str = "CP_and_DP") -> tuple:
    """
    Split a dataframe into metadata and feature data based on the dataset type.

    For CP, only columns containing "CP__" are considered features.
    For DP, only columns containing "DP__" are considered features.
    For CP_and_DP, any column containing "P__" is considered a feature.
    Metadata columns are those that do NOT contain "P__".
    """
    all_cols = df.columns.tolist()
    if dataset == "CP":
        feature_cols = [col for col in all_cols if "CP__" in col]
    elif dataset == "DP":
        feature_cols = [col for col in all_cols if "DP__" in col]
    elif dataset == "CP_and_DP":
        feature_cols = [col for col in all_cols if "P__" in col]
    else:
        raise ValueError(
            "Invalid dataset type. Choose from 'CP', 'DP', or 'CP_and_DP'."
        )

    metadata_cols = [col for col in all_cols if "P__" not in col]
    metadata_dataframe = df[metadata_cols]
    feature_data = df[feature_cols].values
    return metadata_dataframe, feature_data


def load_sc_data(path: Path, is_control: bool = False) -> pd.DataFrame:
    """
    Load single-cell data from a CSV file.

    For negative control data (is_control=True):
      - Insert "Mitocheck_Phenotypic_Class" with value "neg_control".
    """
    df = pd.read_csv(path)
    for c in ["Metadata_Object_Outline", "Unnamed: 0"]:
        if c in df.columns:
            df = df.drop(c, axis=1)
    if is_control:
        df["Mitocheck_Phenotypic_Class"] = "neg_control"
    return df


def merge_profiles(training_df: pd.DataFrame, control_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge training and control single-cell data into one DataFrame.
    """
    merged_df = pd.concat([training_df, control_df], axis=0, ignore_index=True)
    merged_df = merged_df.query(
        "Metadata_Gene != 'failed QC' and Mitocheck_Phenotypic_Class != 'OutOfFocus'"
    )
    return merged_df


def process_merged_cp_features(df: pd.DataFrame, dataset: str = "CP") -> pd.DataFrame:
    """
    Process CP features on the merged dataframe by:
      1. Splitting into metadata and CP feature data.
      2. Concatenating them to form a temporary dataframe.
      3. Running pycytominer's feature_select on the CP features.
      4. Removing original CP columns and appending the feature-selected CP features.
    """
    cp_cols = [col for col in df.columns if col.startswith("CP__")]
    meta, features = split_data(df, dataset=dataset)
    cp_data = pd.concat(
        [meta.reset_index(drop=True), pd.DataFrame(features, columns=cp_cols)], axis=1
    )

    fs_df = feature_select(cp_data, features=cp_cols)
    fs_cp_df = fs_df[[col for col in fs_df.columns if col.startswith("CP__")]]
    df_no_cp = df[[col for col in df.columns if not col.startswith("CP__")]]
    return pd.concat(
        [df_no_cp.reset_index(drop=True), fs_cp_df.reset_index(drop=True)], axis=1
    )


def main():
    training_path = Path("inputs/training_data__no_ic.csv.gz").resolve(strict=True)
    neg_control_path = Path(
        "inputs/normalized_data/negative_control_data.csv.gz"
    ).resolve(strict=True)
    output_dir = Path("outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    training_sc_data = load_sc_data(training_path, is_control=False)
    neg_control_sc_data = load_sc_data(neg_control_path, is_control=True)

    print("Control shape:", neg_control_sc_data.shape)
    print("Training shape:", training_sc_data.shape)

    # Sample negative controls before merging and feature selection.
    neg_control_subset = neg_control_sc_data.sample(frac=0.005, random_state=0)
    print("Negative control subset shape:", neg_control_subset.shape)

    merged_sc_data = merge_profiles(training_sc_data, neg_control_subset)
    print("Merged shape:", merged_sc_data.shape)

    merged_sc_data = process_merged_cp_features(merged_sc_data, dataset="CP")
    print(f"After FS, merged shape: {merged_sc_data.shape}")

    merged_save_path = output_dir / "merged_sc_fs.parquet"
    merged_sc_data.to_parquet(merged_save_path, index=False)
    print(f"Saved merged data to {merged_save_path}")


if __name__ == "__main__":
    main()
