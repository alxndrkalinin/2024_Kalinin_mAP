import numpy as np
import pandas as pd


def generate_plate_map(
    n_wells, n_perts, n_controls, rng, negcon_pert_value="negative_control"
):
    """Generate a mapping of wells to perturbations."""
    wells = [f"w{j}" for j in range(n_wells)]
    perturbations = [f"c{k}" for k in range(n_perts)] + [negcon_pert_value] * n_controls
    shuffled_perts = rng.choice(perturbations, size=len(wells), replace=False)
    return dict(zip(wells, shuffled_perts))


def generate_features_from_distribution(num, params, rng, distribution):
    """Generate features for a given distribution."""
    if distribution == "gaussian":
        return rng.normal(loc=params[0], scale=params[1], size=num)
    elif distribution == "cauchy":
        return rng.standard_cauchy(size=num) * params[1] + params[0]
    else:
        raise ValueError(f"Unknown distribution: {distribution}")


def generate_perturbation_features(
    num_gaussian, params, rng, num_cauchy=0, features_differ=None, differ_params=None
):
    """Generate features for a perturbation using provided distribution parameters."""
    features = []

    if features_differ:
        assert (
            differ_params
        ), "Must provide differ_params if features_differ is not None"

        for dist, num_diff in features_differ.items():
            features.append(
                generate_features_from_distribution(
                    num_diff, differ_params[dist], rng, dist
                )
            )
            if dist == "gaussian":
                num_gaussian -= num_diff
            elif dist == "cauchy":
                num_cauchy -= num_diff

    if num_gaussian > 0:
        features.append(
            generate_features_from_distribution(
                num_gaussian, params["gaussian"], rng, "gaussian"
            )
        )
    if num_cauchy > 0:
        features.append(
            generate_features_from_distribution(
                num_cauchy, params["cauchy"], rng, "cauchy"
            )
        )

    return np.concatenate(features)


def generate_features(
    n_feats,
    n_plates,
    n_wells,
    n_perts,
    n_controls,
    feature_proportions,
    control_params,
    seed=0,
    features_differ=None,
    differ_params=None,
):
    rng = np.random.default_rng(seed)

    num_gaussian = int(n_feats * feature_proportions.get("gaussian", 0))
    num_cauchy = int(n_feats * feature_proportions.get("cauchy", 0))

    plate_map = generate_plate_map(n_wells, n_perts, n_controls, rng)

    metadata_records = []
    all_features = []

    for plate_idx in range(n_plates):
        for well, pert in plate_map.items():
            metadata_records.append(
                {
                    "Metadata_Plate": f"p{plate_idx}",
                    "Metadata_Well": well,
                    "Metadata_Perturbation": pert,
                }
            )

            feature_kwargs = {}
            if pert != "negative_control":
                feature_kwargs = {
                    "features_differ": features_differ,
                    "differ_params": differ_params,
                }

            features = generate_perturbation_features(
                num_gaussian,
                control_params,
                rng,
                num_cauchy=num_cauchy,
                **feature_kwargs,
            )
            all_features.append(features)

    return pd.DataFrame(metadata_records), pd.DataFrame(all_features)
