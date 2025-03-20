import numpy as np
import numpy.typing as npt
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score


def cluster_accuracy(y_true, y_pred):
    accuracy_1 = accuracy_score(y_true, y_pred)
    accuracy_2 = accuracy_score(y_true, np.where(y_pred == 0, 1, 0))
    return max(accuracy_1, accuracy_2)


def kmeans(
    pert_features: npt.ArrayLike,
    control_features: npt.ArrayLike,
    init="random",
    n_clusters: int = 2,
    random_seed: int = 42,
) -> float:
    kmeans = KMeans(n_clusters=n_clusters, init=init, random_state=random_seed)
    merged_features = np.concatenate([pert_features, control_features])
    kmeans.fit(merged_features)
    acc = cluster_accuracy(
        np.concatenate(
            [np.ones(pert_features.shape[0]), np.zeros(control_features.shape[0])]
        ),
        kmeans.labels_,
    )

    return acc
