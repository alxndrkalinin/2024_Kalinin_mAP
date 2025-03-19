import pandas as pd
from pycytominer import aggregate

from map_utils.map import calculate_map
from nelisa_utils import get_meta_features


def load_normalized_profiles():
    """
    Load normalized Cell Painting and nELISA profiles,
    fill missing values, add control indices, and assert equal row counts.
    """
    cp_df = pd.read_parquet("outputs/cellpainting_profiles_normalized.parquet")
    ne_df = pd.read_parquet("outputs/nelisa_profiles_normalized.parquet")

    assert cp_df.shape[0] == ne_df.shape[0], "Row counts do not match!"
    print("Cell Painting shape:", cp_df.shape, "nELISA shape:", ne_df.shape)
    return cp_df, ne_df


def compute_activity_map(
    df, pair_config, map_config, reference_index_config, output_file
):
    """
    Compute phenotypic activity map via calculate_map, filter out 'NA',
    save to CSV, and print the resulting shape.
    """
    result = calculate_map(
        df, pair_config, map_config, reference_index_config=reference_index_config
    )
    result = result.query("Metadata_broad_sample != 'NA'")
    result.to_csv(output_file, index=False)
    print(f"Activity map saved to {output_file} with shape {result.shape}")
    return result


def compute_consistency_map(df, pair_config, map_config):
    """
    Aggregate profiles using pycytominer.aggregate, then compute
    phenotypic consistency map via calculate_map.
    """
    _, feature_cols = get_meta_features(df)
    all_df = aggregate(
        df,
        strata=["Metadata_broad_sample", "Metadata_target_list"],
        features=feature_cols,
    )
    all_df["Metadata_target"] = all_df["Metadata_target_list"].str.split("|")
    all_df = all_df.query("Metadata_broad_sample != 'NA'")
    print("Aggregated consistency DataFrame shape:", all_df.shape)
    result = calculate_map(all_df, pair_config, map_config)
    return result


def main():
    cp_df, ne_df = load_normalized_profiles()

    # Phenotypic activity configuration
    activity_pair_config = {
        "pos_sameby": ["Metadata_broad_sample", "Metadata_control_index"],
        "pos_diffby": [],
        "neg_sameby": [],
        "neg_diffby": ["Metadata_broad_sample", "Metadata_control_index"],
    }
    activity_map_config = {
        "null_size": 100000,
        "groupby_columns": ["Metadata_broad_sample"],
    }

    reference_index_config = {
        "condition": "Metadata_control_type == 'negcon'",
        "reference_col": "Metadata_control_index",
    }
    cp_activity_out = "outputs/cp_map_activity_results.csv"
    ne_activity_out = "outputs/ne_map_activity_results.csv"

    print("Computing phenotypic activity map for Cell Painting...")
    cp_map_results = compute_activity_map(
        cp_df,
        activity_pair_config,
        activity_map_config,
        reference_index_config,
        cp_activity_out,
    )
    cp_map_results.to_csv(cp_activity_out, index=False)
    print(f"Cell Painting activity map saved to {cp_activity_out}")

    print("Computing phenotypic activity map for nELISA...")
    ne_map_results = compute_activity_map(
        ne_df,
        activity_pair_config,
        activity_map_config,
        reference_index_config,
        ne_activity_out,
    )
    ne_map_results.to_csv(ne_activity_out, index=False)
    print(f"nELISA activity map saved to {ne_activity_out}")

    # Phenotypic consistency configuration
    consistency_pair_config = {
        "pos_sameby": ["Metadata_target"],
        "pos_diffby": ["Metadata_broad_sample"],
        "neg_sameby": [],
        "neg_diffby": ["Metadata_target", "Metadata_broad_sample"],
        "multilabel_col": "Metadata_target",
    }
    consistency_map_config = {
        "null_size": 100000,
        "groupby_columns": ["Metadata_target"],
    }
    cp_consistency_out = "outputs/cp_all_map_consistency_results.csv"
    ne_consistency_out = "outputs/ne_all_map_consistency_results.csv"

    print("Computing phenotypic consistency map for Cell Painting...")
    cp_all_map_results = compute_consistency_map(
        cp_df, consistency_pair_config, consistency_map_config
    )
    cp_all_map_results.to_csv(cp_consistency_out, index=False)
    print(f"Cell Painting consistency map saved to {cp_consistency_out}")

    print("Computing phenotypic consistency map for nELISA...")
    ne_all_map_results = compute_consistency_map(
        ne_df, consistency_pair_config, consistency_map_config
    )
    ne_all_map_results.to_csv(ne_consistency_out, index=False)
    print(f"nELISA consistency map saved to {ne_consistency_out}")


if __name__ == "__main__":
    main()
