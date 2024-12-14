from copairs import map as copairs
import numpy as np
import pandas as pd


def merge_labels(drugs_path, samples_path, labels_path, mapper_path):
    drugs = pd.read_csv(drugs_path, sep="\t", skiprows=9)
    samples = pd.read_csv(samples_path, sep="\t", skiprows=9)
    labels = drugs.merge(samples, on="pert_iname")

    # Unify old and new broad_ids
    labels = labels[["deprecated_broad_id", "broad_id", "pert_iname", "target"]]
    labels = labels.melt(id_vars=["pert_iname", "target"]).dropna()
    labels.rename(columns={"value": "broad_id"}, inplace=True)
    labels.drop(columns="variable", inplace=True)

    # Explode broad_ids
    labels["broad_id"] = labels["broad_id"].str.split("|")
    labels = labels.explode("broad_id")
    # Use ignore the last part of the broad_id
    labels["broad_id"] = labels["broad_id"].str[:-9]
    # Explode target
    labels["target"] = labels["target"].str.split("|")
    labels = labels.explode("target")

    # Export broad_id to pert_iname mapper
    mapper = labels.groupby("broad_id", as_index=False).pert_iname.agg(pd.Series.mode)
    mapper.to_csv(mapper_path, index=False)

    # group targets back after deleting dups
    labels = labels.drop_duplicates(subset=["broad_id", "target"])
    labels = labels.groupby("broad_id").agg("|".join)

    # leave samples with single target only
    # labels = labels[~labels['target'].str.contains(r'\|')]
    labels.to_csv(labels_path)


def label_profiles(profiles_path, labels_path, labeled_profiles_path):
    profiles = pd.read_csv(profiles_path, low_memory=False)
    labels = pd.read_csv(labels_path)

    target_map = labels.set_index("broad_id")["target"]

    # ignore the last part of the broad_id
    profiles["Metadata_broad_sample"] = profiles["Metadata_broad_sample"].str[:-9]
    targets = profiles["Metadata_broad_sample"].map(target_map)
    profiles["Metadata_target"] = targets
    profiles = profiles.dropna(subset="Metadata_target")
    # Select a single dose
    profiles = profiles.query("Metadata_dose_recode==6")
    profiles.to_csv(labeled_profiles_path, index=False)


def _split_profile(profiles) -> tuple[pd.DataFrame, np.ndarray]:
    meta_cols = [c for c in profiles if c.startswith("Meta")]
    feat_cols = [c for c in profiles if c not in meta_cols]
    meta = profiles[meta_cols]
    feats = profiles[feat_cols].values
    return meta, feats


def split_consensus_profile(profiles_path) -> tuple[pd.DataFrame, np.ndarray]:
    profiles = pd.read_csv(profiles_path, low_memory=False)
    # Select a single dose
    profiles = profiles.query("Metadata_dose_recode==6")
    return _split_profile(profiles)


def split_profile(profiles_path) -> tuple[pd.DataFrame, np.ndarray]:
    profiles = pd.read_csv(profiles_path, low_memory=False)
    # Select a single dose
    dose_6_plates = profiles.query(  # noqa: F841
        "Metadata_dose_recode==6"
    ).Metadata_Plate.drop_duplicates()
    query = (
        "Metadata_dose_recode==6 or "
        '(Metadata_broad_sample=="DMSO" and Metadata_Plate in @dose_6_plates)'
    )
    profiles = profiles.query(query)
    return _split_profile(profiles)


def _unique_negcons(meta: pd.DataFrame) -> None:
    """
    Hack to avoid mAP computation for negcons. Assign a unique id for every
    negcon so that no pairs are found for such samples.
    """
    negcon_ix = meta["Metadata_broad_sample"] == "DMSO"
    n_negcon = negcon_ix.sum()
    negcon_ids = [f"DMSO_{i}" for i in range(n_negcon)]
    meta.loc[negcon_ix, "Metadata_broad_sample"] = negcon_ids


def ap_technical(profiles_path, ap_path):
    meta, feats = split_profile(profiles_path)
    _unique_negcons(meta)

    result = copairs.average_precision(
        meta,
        feats,
        pos_sameby=["Metadata_broad_sample"],
        # 'Metadata_Well' is not used because replicates are in the same position
        pos_diffby=[],
        neg_sameby=["Metadata_Plate"],
        neg_diffby=["Metadata_pert_type", "Metadata_broad_sample"],
    )
    result.to_csv(ap_path, index=False)


def ap_biological(profiles_path, ap_path):
    meta, feats = split_consensus_profile(profiles_path)

    # Filter out DMSO
    ix = meta["Metadata_broad_sample"] != "DMSO"
    meta = meta.loc[ix].reset_index(drop=True).copy()
    feats = feats[ix]

    # split target labels
    meta["Metadata_target"] = meta["Metadata_target"].str.split("|")

    result = copairs.multilabel.average_precision(
        meta,
        feats,
        pos_sameby=["Metadata_target"],
        pos_diffby=["Metadata_broad_sample"],
        neg_sameby=[],
        neg_diffby=["Metadata_target"],
        multilabel_col="Metadata_target",
    )
    result.to_csv(ap_path, index=False)


def mean_average_precision(ap_path, map_path, sameby, pvalue_threshold):
    ap_scores = pd.read_csv(ap_path)
    agg_result = copairs.mean_average_precision(
        ap_scores, sameby, null_size=10000, threshold=pvalue_threshold, seed=0
    )
    agg_result.to_csv(map_path, index=False)
