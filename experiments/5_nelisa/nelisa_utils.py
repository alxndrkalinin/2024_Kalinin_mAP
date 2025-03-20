from functools import partial

from pycytominer import feature_select, normalize


def get_meta_features(df):
    meta_cols = df.filter(regex="Metadata_").columns.tolist()
    feat_cols = df.filter(regex="^(?!Metadata)").columns.tolist()
    return meta_cols, feat_cols


def normalize_select(
    df, norm_group_col="Metadata_Plate", select_features=None, inplace=True
):
    """
    Normalize and select features
    """
    df = df.copy() if not inplace else df

    meta_cols, feat_cols = get_meta_features(df)
    norm_func = partial(
        normalize, method="mad_robustize", features=feat_cols, meta_features=meta_cols
    )
    df = df.groupby(norm_group_col).apply(norm_func).reset_index(drop=True)
    if select_features is not None:
        df = feature_select(df, operation=select_features, na_cutoff=0)
    return df


def unify_pert_ids(df1, df2, pert_col="Metadata_broad_sample"):
    """
    Unify perturbation ids between two dataframes
    """
    perts1 = set(df1[pert_col].unique())
    perts2 = set(df2[pert_col].unique())
    diff_perts = list((perts1 | perts2) - (perts1 & perts2))  # noqa: F841
    df1 = df1.query(f"{pert_col} not in @diff_perts")
    df2 = df2.query(f"{pert_col} not in @diff_perts")

    counts_diff = (
        df1.groupby("Metadata_broad_sample").size()
        != df2.groupby("Metadata_broad_sample").size()
    )
    diff_count_perts = counts_diff.index[counts_diff.values]  # noqa: F841
    df1 = df1.query(f"{pert_col} not in @diff_count_perts").reset_index(drop=True)
    df2 = df2.query(f"{pert_col} not in @diff_count_perts").reset_index(drop=True)

    return df1, df2
