import numpy as np

from copairs import map


def assign_reference_index(
    df, condition, reference_col="Metadata_Reference_Index", default_value=-1
):
    """Assign reference index to a column based on a condition."""
    df[reference_col] = default_value
    df.loc[condition, reference_col] = df.loc[condition].index
    return df


def calculate_map(df, pair_config, map_config):
    """Calculate mean average precision (mAP) for a given DataFrame."""
    metadata = df.filter(regex="Metadata_")
    features = df.filter(regex="^(?!Metadata)").values

    pos_sameby = pair_config.get("pos_sameby", {"all": [], "any": []})
    pos_diffby = pair_config.get("pos_diffby", {"all": [], "any": []})
    neg_sameby = pair_config.get("neg_sameby", {"all": [], "any": []})
    neg_diffby = pair_config.get("neg_diffby", {"all": [], "any": []})
    multilabel_col = pair_config.get("multilabel_col", None)

    if multilabel_col is not None:
        ap_scores = map.multilabel.average_precision(
            metadata,
            features,
            pos_sameby,
            pos_diffby,
            neg_sameby,
            neg_diffby,
            multilabel_col,
        )
    else:
        ap_scores = map.average_precision(
            metadata, features, pos_sameby, pos_diffby, neg_sameby, neg_diffby
        )
    ap_scores.dropna(subset=["average_precision"], inplace=True)

    null_size = map_config.get("null_size", 10000)
    p_value_threshold = map_config.get("threshold", 0.05)
    random_seed = map_config.get("random_seed", 0)
    groupby_columns = map_config.get("groupby_columns", pos_sameby)
    map_scores = map.mean_average_precision(
        ap_scores,
        groupby_columns,
        null_size=null_size,
        threshold=p_value_threshold,
        seed=random_seed,
    )

    map_scores["-log10(mAP p-value)"] = -np.log10(
        map_scores["corrected_p_value"].clip(lower=np.finfo(float).eps)
    )
    map_scores.rename(
        columns={"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        inplace=True,
    )
    return map_scores
