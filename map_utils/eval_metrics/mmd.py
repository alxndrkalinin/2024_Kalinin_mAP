from typing import Optional

import numpy as np
import numpy.typing as npt
from scipy.spatial.distance import squareform, pdist, cdist


def mmd(
    pert_features: npt.ArrayLike,
    control_features: npt.ArrayLike,
    kernel: str = "rbf",
    sigma: Optional[float] = None,
    nb_permutations: int = 1000,
    random_seed: int = 42,
):
    mmd = quadratic_time_mmd(
        x=pert_features, y=control_features, kernel=kernel, sigma=sigma
    )

    mmd_sim = np.zeros(nb_permutations)
    rng = np.random.default_rng(random_seed)

    merged_features = np.concatenate([pert_features, control_features])
    pert_mask = np.concatenate(
        [np.ones(pert_features.shape[0]), np.zeros(control_features.shape[0])]
    )

    for i in range(nb_permutations):
        pert_mask_perm = rng.permutation(pert_mask).astype(bool)
        mmd_sim[i] = quadratic_time_mmd(
            x=merged_features[pert_mask_perm],
            y=merged_features[np.logical_not(pert_mask_perm)],
            kernel="rbf",
            sigma=None,
        )

    return (np.sum(mmd_sim >= mmd) + 1) / (nb_permutations + 1)


def sq_distances(x, y=None):
    """
    If Y=None, then this computes the distance between X and itself
    """
    assert x.ndim == 2

    if y is None:
        sq_dists = squareform(pdist(x, "sqeuclidean"))
    else:
        assert y.ndim == 2
        assert x.shape[1] == y.shape[1]
        sq_dists = cdist(x, y, "sqeuclidean")

    return sq_dists


def gauss_kernel(x, y=None, sigma=1.0):
    """
    Computes the standard Gaussian kernel k(x,y)=exp(- ||x-y||**2 / (2 * sigma**2))

    X - 2d array, samples on left hand side
    Y - 2d array, samples on right hand side, can be None in which case they are replaced by X

    returns: kernel matrix
    """
    sq_dists = sq_distances(x, y)
    K = np.exp(-sq_dists / (2 * sigma**2))
    return K


def linear_kernel(x, y):
    return np.dot(x, y.T)


def gaussian_kernel_median_heuristic(Z):
    sq_dists = sq_distances(Z)
    np.fill_diagonal(sq_dists, np.nan)
    sq_dists = np.ravel(sq_dists)
    sq_dists = sq_dists[~np.isnan(sq_dists)]
    median_dist = np.median(np.sqrt(sq_dists))
    return np.sqrt(median_dist / 2.0)


def quadratic_time_mmd(x, y, kernel="rbf", sigma=None):
    assert x.ndim == y.ndim == 2

    if kernel == "rbf":
        if sigma is None:
            sigma = gaussian_kernel_median_heuristic(np.vstack([x, y]))

        def kernel_fn(a, b):
            return gauss_kernel(a, b, sigma)
    elif kernel == "linear":
        kernel_fn = linear_kernel
    else:
        raise ValueError("Unknown kernel")

    k_xx = kernel_fn(x, x)
    k_xy = kernel_fn(x, y)
    k_yy = kernel_fn(y, y)

    n = len(k_xx)
    m = len(k_yy)

    np.fill_diagonal(k_xx, 0)
    np.fill_diagonal(k_yy, 0)
    mmd = (
        np.sum(k_xx) / (n * (n - 1))
        + np.sum(k_yy) / (m * (m - 1))
        - 2 * np.sum(k_xy) / (n * m)
    )
    return mmd
