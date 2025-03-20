"""
Mostly borrowed from https://github.com/jump-cellpainting/2024_Chandrasekaran_Morphmap/blob/main/03.retrieve-annotations/utils.py

> Chandrasekaran, S. N. et al. Morphological map of under- and over-expression of genes in human cells. bioRxiv 2024.12.02.624527 (2024)
"""

import numpy as np
import pandas as pd


def get_metacols(df):
    """return a list of metadata columns"""
    return [c for c in df.columns if c.startswith("Metadata_")]


def get_featurecols(df):
    """returna  list of featuredata columns"""
    return [c for c in df.columns if not c.startswith("Metadata")]


def get_metadata(df):
    """return dataframe of just metadata columns"""
    return df[get_metacols(df)]


def get_featuredata(df):
    """return dataframe of just featuredata columns"""
    return df[get_featurecols(df)]


def remove_low_infection_efficiency_wells(df):
    """return dataframe of non-low infection efficiency wells"""

    reagents_to_remove_df = pd.read_csv(
        "inputs/orf-reagents-low-infection-efficiency-and-outliers.csv.gz",
        usecols=["Metadata_plate_map_name", "Metadata_Well", "Metadata_broad_sample"],
    )

    df = (
        df.merge(
            reagents_to_remove_df,
            on=["Metadata_plate_map_name", "Metadata_Well", "Metadata_broad_sample"],
            how="left",
            indicator=True,
        )
        .query('_merge=="left_only"')
        .drop(columns=["_merge"])
    )

    return df


def remove_empty_wells(df):
    """return dataframe of non-empty wells"""
    df = df.query("Metadata_pert_type!='empty'").reset_index(drop=True)
    df = df.dropna(subset=["Metadata_broad_sample"]).reset_index(drop=True)
    df = df.dropna(subset=["Metadata_Symbol"]).reset_index(drop=True)
    return df


def consensus(profiles_df, group_by_feature):
    """
    Computes the median consensus profiles.
    Parameters:
    -----------
    profiles_df: pandas.DataFrame
        dataframe of profiles
    group_by_feature: str
        Name of the column
    Returns:
    -------
    pandas.DataFrame
    """

    metadata_df = get_metadata(profiles_df).drop_duplicates(subset=[group_by_feature])

    feature_cols = [group_by_feature] + get_featurecols(profiles_df)
    profiles_df = (
        profiles_df[feature_cols].groupby([group_by_feature]).median().reset_index()
    )

    profiles_df = metadata_df.merge(profiles_df, on=group_by_feature)

    return profiles_df


def remove_nan_features(df):
    """
    Remove NaN features.
    Parameters
    ----------
    df: pandas.DataFrame
        DataFrame of profiles.
    Returns
    -------
    df: pandas.DataFrame
        DataFrame of profiles after dropping NaN features.
    """

    r, c = np.where(df.isna())
    features_to_remove = [
        _ for _ in list(df.columns[list(set(c))]) if not _.startswith("Metadata_")
    ]
    print(f"Removed nan features: {features_to_remove}")
    df.drop(features_to_remove, axis=1, inplace=True)
    return df
