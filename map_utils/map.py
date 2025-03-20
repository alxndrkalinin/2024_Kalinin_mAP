import numpy as np
from statsmodels.stats.multitest import multipletests

from copairs import map, matching
from copairs.map.average_precision import p_values


def calculate_map(
    df, pair_config, map_config, reference_index_config=None, return_ap_scores=False
):
    """Calculate mean average precision (mAP) for a given DataFrame."""
    if reference_index_config is not None:
        df = matching.assign_reference_index(df, **reference_index_config)
    metadata = df.filter(regex="Metadata_")
    features = df.filter(regex="^(?!Metadata)").values

    pos_sameby = pair_config.get("pos_sameby", [])
    pos_diffby = pair_config.get("pos_diffby", [])
    neg_sameby = pair_config.get("neg_sameby", [])
    neg_diffby = pair_config.get("neg_diffby", [])
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
    max_workers = map_config.get("max_workers", 32)
    p_value_threshold = map_config.get("threshold", 0.05)
    random_seed = map_config.get("random_seed", 0)
    groupby_columns = map_config.get("groupby_columns", pos_sameby)
    map_scores = map.mean_average_precision(
        ap_scores,
        groupby_columns,
        null_size=null_size,
        threshold=p_value_threshold,
        seed=random_seed,
        max_workers=max_workers,
    )

    map_scores["-log10(mAP p-value)"] = -np.log10(
        map_scores["corrected_p_value"].clip(lower=np.finfo(float).eps)
    )
    map_scores.rename(
        columns={"mean_average_precision": "mAP", "below_corrected_p": "p < 0.05"},
        inplace=True,
    )

    if return_ap_scores:
        ap_p_values = p_values(ap_scores, null_size=null_size, seed=random_seed)
        _, ap_p_values, _, _ = multipletests(ap_p_values, method="fdr_bh")
        ap_scores["p_value"] = ap_p_values
        ap_scores["p < 0.05"] = ap_p_values < p_value_threshold
        ap_scores["-log10(AP p-value)"] = -np.log10(
            ap_p_values.clip(min=np.finfo(float).eps)
        )
        return map_scores, ap_scores
    else:
        return map_scores
