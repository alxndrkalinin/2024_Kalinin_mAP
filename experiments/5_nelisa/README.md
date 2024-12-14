## 5. nELISA dataset

This experiment demostrated application of mAP to single-cell proteomic profiling data.

### Dataset description

We used the dataset containing proteomic profiles from a 191-plex nELISA, a high-throughput, high-plex assay designed for quantitative profiling of the secretome, which was performed in A549 cells across 306 well-characterized compound perturbations from the Broad Institute’s drug repurposing library. This dataset also included matching CellProfiler morphological profiles from Cell Painting images of the same physical samples whose supernatants were nELISA-profiled. Profiles were normalized per-plate by subtracting medians from feature values and dividing them by median absolute deviation (“Robust MAD”). Feature selection included variance thresholding to remove features with minimal variation across the dataset, removing highly correlated features, and removing features with missing values or outliers.

Citation:
> Dagher, M. et al. nELISA: A high-throughput, high-plex platform enables quantitative profiling of the secretome. bioRxiv (2023) doi:10.1101/2023.04.17.535914.

## How to run

nELISA and Cell Painting profiles are provided within `inputs` directory.

0. [Preprocess profiles](./0_preprocess_profiles.ipynb) by performing robust normalization and feature selection.
1. [Calculate mAP](./1_calculate_map.ipynb) to assess phenotypic activity and consistency.
2. [Plot results](./2_plot_map.ipynb).
