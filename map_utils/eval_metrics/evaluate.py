"""Calculate evaluation metrics from profiling experiments.

Loosely based on now-deprecated cytominer-eval package:
https://github.com/cytomining/cytominer-eval.
"""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from itertools import chain
from typing import List, Optional, Tuple, Union

import networkx as nx
import numpy.typing as npt
import pandas as pd
from tqdm.auto import tqdm

from copairs.map.filter import evaluate_and_filter, extract_filters, flatten_str_list
from copairs.matching import Matcher, MatcherMultilabel, dict_to_dframe
from map_utils.eval_metrics import kmeans, mean_ap, mmd, mp_value


def calculate_accuracy(metrics_results: dict) -> dict:
    """Calculate accuracy from metric results."""
    acc = {}
    for metric, results in metrics_results.items():
        if metric == "kmeans":
            acc[metric] = (results["kmeans"] == 1).mean()
        else:
            acc[metric] = (results[metric] < 0.05).mean()
    return acc


def get_copairs(
    meta: pd.DataFrame,
    pos_sameby: List[str],
    pos_diffby: List[str],
    neg_sameby: List[str],
    neg_diffby: List[str],
    filters: Optional[List[str]] = None,
    multilabel_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Get positive and negative pairs."""
    meta = meta.reset_index(drop=True).copy()

    # generic filters that do not affect matching
    if filters is not None:
        meta, columns = evaluate_and_filter(meta, filters)

    columns = flatten_str_list(pos_sameby, pos_diffby, neg_sameby, neg_diffby)
    if all(c in meta.columns for c in columns):
        matcher = create_matcher(meta, columns, multilabel_col)

        pos_pairs = matcher.get_all_pairs(sameby=pos_sameby, diffby=pos_diffby)
        pos_pairs = dict_to_dframe(pos_pairs, flatten_str_list(pos_sameby))

        neg_pairs = matcher.get_all_pairs(sameby=neg_sameby, diffby=neg_diffby)
        neg_pairs = set(chain.from_iterable(neg_pairs.values()))
        neg_pairs = pd.DataFrame(neg_pairs, columns=["ix1", "ix2"])

    else:
        # filtering the whole meta df may remove raws needed for negative pairs
        #  so we get positive and negative pairs separately
        pos_columns = flatten_str_list(pos_sameby, pos_diffby)
        matcher = create_matcher(meta, pos_columns, multilabel_col)
        pos_pairs = matcher.get_all_pairs(sameby=pos_sameby, diffby=pos_diffby)
        _, pos_sameby = extract_filters(flatten_str_list(pos_sameby), meta.columns)
        pos_pairs = dict_to_dframe(pos_pairs, pos_sameby)

        neg_columns = flatten_str_list(neg_sameby, neg_diffby)
        matcher = create_matcher(meta, neg_columns, multilabel_col)
        neg_pairs = matcher.get_all_pairs(sameby=neg_sameby, diffby=neg_diffby)
        neg_pairs = set(chain.from_iterable(neg_pairs.values()))
        neg_pairs = pd.DataFrame(neg_pairs, columns=["ix1", "ix2"])
        print(
            f"Pos pairs size: {pos_pairs.shape[0]}, Neg pairs size: {neg_pairs.shape[0]}"
        )

    return pos_pairs, neg_pairs


def create_matcher(
    meta: pd.DataFrame, columns: list, multilabel_col=None
) -> Union[Matcher, MatcherMultilabel]:
    """Create a matcher object."""
    meta, columns = evaluate_and_filter(meta, columns)
    if multilabel_col is None:
        return Matcher(meta, columns, seed=0)
    else:
        return MatcherMultilabel(meta, columns, multilabel_col, seed=0)


def get_pos_neg_pairs(metadata: pd.DataFrame, pair_grouping: dict):
    """Get positive and negative pairs."""
    pos_sameby = pair_grouping.get("pos_sameby", [])
    pos_diffby = pair_grouping.get("pos_diffby", [])
    neg_sameby = pair_grouping.get("neg_sameby", [])
    neg_diffby = pair_grouping.get("neg_diffby", [])
    multilabel_col = pair_grouping.get("multilabel_col", None)

    pair_grouping = flatten_str_list(pos_sameby, pos_diffby, neg_sameby, neg_diffby)

    pos_pairs, neg_pairs = get_copairs(
        meta=metadata,
        pos_sameby=pos_sameby,
        pos_diffby=pos_diffby,
        neg_sameby=neg_sameby,
        neg_diffby=neg_diffby,
        multilabel_col=multilabel_col,
    )

    return pos_pairs, neg_pairs


def map_neg_pos_cliques(pos_cliques: npt.ArrayLike, neg_cliques: npt.ArrayLike) -> dict:
    """Map negative cliques to positive cliques."""
    # create a reverse lookup dictionary for neg_cliques
    neg_lookup = defaultdict(list)
    for nc in neg_cliques:
        for elem in nc:
            neg_lookup[elem].extend(nc)

    # flatten the lists in neg_lookup and remove duplicates
    for key in neg_lookup:
        neg_lookup[key] = list(set(neg_lookup[key]))

    clique_dict = {}
    for index, pos_clique in enumerate(pos_cliques):
        # collect all unique neg_clique elements for each pos_clique
        neg_elements = set()
        for elem in pos_clique:
            neg_elements.update(neg_lookup.get(elem, []))
        # remove any element from neg_elements that is in pos_clique to avoid self-reference
        clique_dict[index] = sorted(neg_elements.difference(pos_clique))

    return clique_dict


def process_group(
    group, grouped_profiles, grouped_controls, features, metric, metric_kwargs
) -> Tuple[int, float]:
    """Process a group of profiles."""
    group_id, group_df = group
    if group_id in grouped_controls:
        control_df = grouped_profiles.loc[grouped_controls[group_id], :]
    else:
        return group_id, pd.NA

    metric_fn = get_metric_fn(metric)
    metric_value = metric_fn(
        group_df[features].values,
        control_df[features].values,
        **metric_kwargs,
    )
    return group_id, metric_value


def process_groups(
    grouped_profiles: pd.DataFrame,
    features: list,
    metric: str,
    group_col: str = "group_id",
    metric_kwargs: Optional[dict] = None,
) -> pd.DataFrame:
    """Process profiles by groups."""
    profiles, controls = grouped_profiles
    partial_process_group = partial(
        process_group,
        grouped_profiles=profiles,
        grouped_controls=controls,
        features=features,
        metric=metric,
        metric_kwargs=metric_kwargs,
    )

    with ThreadPoolExecutor(max_workers=1) as executor:  # debug
        # with ThreadPoolExecutor() as executor:
        profile_groups = list(
            tqdm(profiles.groupby(group_col), desc=f"Calculating {metric}")
        )
        results = list(
            tqdm(
                executor.map(partial_process_group, profile_groups),
                total=len(profile_groups),
                desc="Processing groups",
            )
        )

    results = pd.DataFrame(results, columns=[group_col, f"{metric}"])
    return results.reset_index(drop=True)


def group_profiles(
    profiles: pd.DataFrame, metadata_features: list, metadata_groups: dict
) -> Tuple[pd.DataFrame, dict]:
    """Group profiles in cliques based on positive and negative pairs."""
    metadata = profiles.loc[:, metadata_features]
    pos_pairs, neg_pairs = get_pos_neg_pairs(metadata, metadata_groups)

    pos_graph = nx.Graph()
    pos_graph.add_edges_from(pos_pairs[["ix1", "ix2"]].values.tolist())
    pos_cliques = list(nx.find_cliques(pos_graph))
    assert len(set.intersection(*map(set, pos_cliques))) == 0, (
        "Error! Overlapping positive cliques found."
    )

    neg_graph = nx.Graph()
    neg_graph.add_edges_from(neg_pairs[["ix1", "ix2"]].values.tolist())
    neg_cliques = list(nx.find_cliques(neg_graph))
    pos_to_neg_map = map_neg_pos_cliques(pos_cliques, neg_cliques)

    grouped_profiles = profiles.copy()
    grouped_profiles["group_id"] = pd.NA
    for idx, clique in enumerate(pos_cliques):
        grouped_profiles.loc[clique, "group_id"] = idx
        grouped_profiles.loc[pos_to_neg_map[idx], "group_id"] = -1

    return grouped_profiles.dropna(axis=0, subset=["group_id"]), pos_to_neg_map


def get_metric_fn(operation_name: str) -> callable:
    """Get the function for a given operation name."""
    operation_mapping = {
        "mp_value": mp_value,
        "mean_ap": mean_ap,
        "kmeans": kmeans,
        "mmd": mmd,
    }
    try:
        return operation_mapping[operation_name]
    except KeyError:
        raise ValueError(f"Invalid operation name: {operation_name}")


def evaluate_metric(
    profiles: pd.DataFrame,
    features: List[str],
    metadata_features: List[str],
    metadata_groups: Union[List[str], dict],
    metric: str,
    metric_kwargs: dict,
    grouped_profiles: Optional[Tuple[pd.DataFrame]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate a metric on grouped profiles."""

    if grouped_profiles is None:
        grouped_profiles = group_profiles(profiles, metadata_features, metadata_groups)

    print(f"\nCalculating metric: {metric}")
    metric_result = process_groups(
        grouped_profiles, features, metric, metric_kwargs=metric_kwargs
    )

    pos_meta_cols = flatten_str_list(
        metadata_groups.get("pos_sameby", []), metadata_groups.get("pos_diffby", [])
    )

    metadata = grouped_profiles[0].loc[:, ~grouped_profiles[0].columns.isin(features)]
    _, pos_meta_cols = extract_filters(pos_meta_cols, metadata.columns)

    metadata = metadata[pos_meta_cols + ["group_id"]]
    metric_result = (
        metadata.merge(metric_result, on="group_id", how="left")
        .dropna(subset=[metric])
        .drop(columns="group_id")
        .drop_duplicates()
        .reset_index(drop=True)
    )
    metric_result = metric_result[pos_meta_cols + [metric]]

    return metric_result, grouped_profiles


def evaluate_metrics(
    profiles: pd.DataFrame,
    features: List[str],
    metadata_features: List[str],
    metadata_groups: Union[List[str], dict],
    metrics_config: dict,
    grouped_profiles: Optional[pd.DataFrame] = None,
) -> dict:
    """Group profiles and evaluate metrics."""
    metric_results = {}
    for metric, metric_kwargs in metrics_config.items():
        metric_result = evaluate_metric(
            profiles=profiles,
            features=features,
            metadata_features=metadata_features,
            metadata_groups=metadata_groups,
            metric=metric,
            metric_kwargs=metric_kwargs,
            grouped_profiles=grouped_profiles,
        )
        metric_results[metric] = metric_result[0]
        grouped_profiles = metric_result[1]

    return metric_results
