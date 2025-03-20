import pandas as pd
from pathlib import Path
from nelisa_utils import normalize_select, unify_pert_ids


def load_and_normalize(norm_confis: dict) -> dict:
    """
    Load profiles from parquet files and normalize them according to the provided configurations.

    Parameters
    ----------
    norm_confis : dict
        Dictionary with normalization configurations for each dataset.

    Returns
    -------
    dict
        Dictionary of normalized DataFrames keyed by dataset name.
    """
    normed_dfs = {}
    for norm_name, norm_conf in norm_confis.items():
        path = norm_conf.pop("path")
        df = pd.read_parquet(path)
        print(f"Normalizing {norm_name} profiles of shape {df.shape}...")
        normed_df = normalize_select(df, **norm_conf)
        print(f"Normalized {norm_name} profiles to shape {normed_df.shape}")
        normed_df.fillna({"Metadata_broad_sample": "NA"}, inplace=True)
        normed_dfs[norm_name] = normed_df
    return normed_dfs


def unify_profiles(normed_dfs: dict, pert_col: str = "Metadata_broad_sample") -> dict:
    """
    Unify perturbation IDs across datasets.

    Parameters
    ----------
    normed_dfs : dict
        Dictionary of normalized DataFrames.
    pert_col : str, optional
        Column name to use for unification, by default "Metadata_broad_sample".

    Returns
    -------
    dict
        Dictionary with unified DataFrames.
    """
    unified_dfs = unify_pert_ids(*list(normed_dfs.values()), pert_col=pert_col)
    # Check that all unified DataFrames have the same number of rows
    assert all(len(df) == len(unified_dfs[0]) for df in unified_dfs), (
        "Unified DataFrames have mismatched lengths."
    )

    for idx, norm_name in enumerate(normed_dfs.keys()):
        print(f"Unified {norm_name} profiles to shape {unified_dfs[idx].shape}")
        normed_dfs[norm_name] = unified_dfs[idx]
    return normed_dfs


def save_profiles(normed_dfs: dict, output_dir: Path = Path("outputs")) -> None:
    """
    Save normalized profiles to parquet files.

    Parameters
    ----------
    normed_dfs : dict
        Dictionary of normalized DataFrames.
    output_dir : Path, optional
        Directory to save the files, by default "outputs".
    """
    output_dir.mkdir(exist_ok=True, parents=True)
    for norm_name, normed_df in normed_dfs.items():
        out_file = output_dir / f"{norm_name}_profiles_normalized.parquet"
        normed_df.to_parquet(out_file)
        print(f"Saved normalized {norm_name} profiles to {out_file}")


def main():
    # Configuration for feature selection in Cell Painting profiles
    cp_feature_select_ops = [
        "variance_threshold",
        "correlation_threshold",
        "drop_na_columns",
        "blocklist",
        "drop_outliers",
    ]

    norm_confis = {
        "cellpainting": {
            "path": "inputs/cellpainting_profiles.parquet",
            "norm_group_col": "Metadata_Plate",
            "select_features": cp_feature_select_ops,
        },
        "nelisa": {
            "path": "inputs/nelisa_profiles.parquet",
            "norm_group_col": "Metadata_nelisa_plate_id",
            "select_features": None,
        },
    }
    normed_dfs = load_and_normalize(norm_confis)
    normed_dfs = unify_profiles(normed_dfs, pert_col="Metadata_broad_sample")
    save_profiles(normed_dfs)


if __name__ == "__main__":
    main()
