## 7. Mitocheck datasets

This experiment demostrated application of mAP to single-cell image-based morphological profiles.

### Dataset description

We used the previously published Mitocheck dataset containing images of GFP-tagged nuclei of HeLa cells perturbed with small interfering RNA (siRNAs) to silence approximately 21,000 protein-coding genes. Within the dataset, approximately 3,000 cell images were manually labeled into one of 15 morphological phenotype classes. Recently, these images were re-analyzed with a more comprehensive image analysis pipeline, which included illumination correction using PyBasic, segmentation using CellPose, and single-cell feature extraction using CellProfiler and DeepProfiler. Extracted profiles were standardized by removing the mean and scaling to the unit variance of negative control cells. We performed feature selection for both CellProfiler- and DeepProfiler-derived profiles by variance thresholding to remove features with minimal variation across the dataset, removing highly correlated features, removing features with missing values or outliers, and removing blocklisted features—all using pycytominer default parameters.

Imaging data reference:
> Neumann, B. et al. Phenotypic profiling of the human genome by time-lapse microscopy reveals cell division genes. Nature 464, 721–727 (2010).

Profiles reference:
> Tomkinson, J., Kern, R., Mattson, C. & Way, G. P. Toward generalizable phenotype prediction from single-cell morphology representations. BMC Methods 1, 1–14 (2024).

### How to run

0. Download profiling data by running `python 0_download_data.py`.
1. Preprocess data with `python 0_prepare_data.py`.
2. Calculate mAP with `python 2_calculate_map.py`.
3. [Plot results](./3_plot_map.ipynb).
