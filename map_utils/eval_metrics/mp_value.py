"""Functions to calculate multidimensional perturbation values (mp-value)

mp-value describes the distance, in dimensionality-reduced space, between a perturbation
and a control [1]_.

References
----------

.. [1] Hutz, J. et al. "The Multidimensional Perturbation Value: A Single Metric to
   Measure Similarity and Activity of Treatments in High-Throughput Multidimensional
   Screens" Journal of Biomolecular Screening, Volume: 18 issue: 4, page(s): 367-377.
   doi: 10.1177/1087057112469257
"""

from typing import Union

import numpy as np
import pandas as pd
import numpy.typing as npt
from sklearn.decomposition import PCA
from sklearn.covariance import EmpiricalCovariance


class MahalanobisEstimator:
    """
    Store location and dispersion estimators of the empirical distribution of data
    provided in an array and allow computation of statistical distances.

    Parameters
    ----------
    arr : {pandas.DataFrame, np.ndarray}
        the matrix used to calculate covariance

    Attributes
    ----------
    sigma : np.array
        Fitted covariance matrix of sklearn.covariance.EmpiricalCovariance()

    Methods
    -------
    mahalanobis(X)
        Computes mahalanobis distance between the input array (self.arr) and the X
        array as provided
    """

    def __init__(self, arr: Union[pd.DataFrame, np.ndarray]):
        self.sigma = EmpiricalCovariance().fit(arr)

    def mahalanobis(self, X: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
        """Compute the mahalanobis distance between the empirical distribution described
        by this object and points in an array `X`.

        Parameters
        ----------
        X : {pandas.DataFrame, np.ndarray}
            A samples by features array-like matrix to compute mahalanobis distance
            between self.arr

        Returns
        -------
        numpy.array
            Mahalanobis distance between the input array and the original sigma
        """
        return self.sigma.mahalanobis(X)


def calculate_mahalanobis(
    pert_features: npt.ArrayLike, control_features: npt.ArrayLike
) -> npt.ArrayLike:
    """Given perturbation and control dataframes, calculate mahalanobis distance per
    perturbation

    Usage: Designed to be called within a pandas.DataFrame().groupby().apply(). See
    :py:func:`cytominer_eval.operations.util.calculate_mp_value`.

    Parameters
    ----------
    pert_df : pandas.DataFrame
        A pandas dataframe of replicate perturbations (samples by features)
    control_df : pandas.DataFrame
        A pandas dataframe of control perturbations (samples by features). Must have the
        same feature measurements as pert_df

    Returns
    -------
    float
        The mahalanobis distance between perturbation and control
    """
    assert len(control_features) > 1, "Error! No control perturbations found."

    # Get dispersion and center estimators for the control perturbations
    control_estimators = MahalanobisEstimator(control_features)

    # Distance between mean of perturbation and control
    maha = control_estimators.mahalanobis(
        np.array(np.mean(pert_features, 0)).reshape(1, -1)
    )[0]
    return maha


def mp_value(
    pert_features: npt.ArrayLike,
    control_features: npt.ArrayLike,
    rescale_pca: bool = True,
    nb_permutations: int = 1000,
    random_seed: int = 42,
) -> npt.ArrayLike:
    """Given perturbation and control dataframes, calculate mp-value per perturbation

    Usage: Designed to be called within a pandas.DataFrame().groupby().apply(). See
    :py:func:`cytominer_eval.operations.mp_value.mp_value`.

    Parameters
    ----------
    pert_df : pandas.DataFrame
        A pandas dataframe of replicate perturbations (samples by features)
    control_df : pandas.DataFrame
        A pandas dataframe of control perturbations (samples by features). Must have the
        same feature measurements as pert_df
    params : {dict}, optional
        the parameters to use when calculating mp value. See
        :py:func:`cytominer_eval.operations.util.default_mp_value_parameters`.

    Returns
    -------
    float
        The mp value for the given perturbation

    """
    assert len(control_features) > 1, "Error! No control perturbations found."

    merged_features = np.concatenate([pert_features, control_features])

    # We reduce the dimensionality with PCA
    # so that 90% of the variance is conserved
    pca = PCA(n_components=0.9, svd_solver="full", random_state=random_seed)

    try:
        pca_array = pca.fit_transform(merged_features)
    except np.linalg.LinAlgError as err:
        if "SVD did not converge" in str(err):
            print(
                "SVD did not converge: check that merged dataframe does not contain duplicate rows or columns."
            )
        raise err

    # We scale columns by the variance explained
    if rescale_pca:
        pca_array = pca_array * pca.explained_variance_ratio_
    # This seems useless, as the point of using the Mahalanobis
    # distance instead of the Euclidean distance is to be independent
    # of axes scales

    # Distance between mean of perturbation and control
    obs = calculate_mahalanobis(
        pca_array[: pert_features.shape[0]],
        pca_array[-control_features.shape[0] :],
    )

    # Permutation test
    sim = np.zeros(nb_permutations)
    pert_mask = np.concatenate(
        [np.ones(pert_features.shape[0]), np.zeros(control_features.shape[0])]
    )

    rng = np.random.default_rng(random_seed)
    for i in range(nb_permutations):
        pert_mask_perm = rng.permutation(pert_mask).astype(bool)
        pert_perm = pca_array[pert_mask_perm]
        control_perm = pca_array[np.logical_not(pert_mask_perm)]
        sim[i] = calculate_mahalanobis(pert_perm, control_perm)

    return (np.sum(sim >= obs) + 1) / (nb_permutations + 1)
