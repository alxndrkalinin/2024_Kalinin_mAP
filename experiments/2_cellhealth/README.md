## 2. Cell Health

This experiment demostrated versatility of the mAP framework through its application to different tasks on “Cell Health” dataset, evaluating the effects of selected preprocessing methods and experimental designs.

### Dataset description

We used our previously published “Cell Health” dataset of Cell Painting images of CRISPR-Cas9 knockout perturbations of 59 genes, targeted by 119 guides in three different cell lines (A549, ES2, and HCC44). Morphological profiles were previously extracted from images using CellProfiler47 and median-aggregated on the well level. We used a subset of 100 guides that had exactly six replicates (two replicates in three different plates) in each cell line. We performed two types of profile normalization followed by feature selection using pycytominer. The first preprocessing method included data standardization by subtracting means from feature values and dividing them by variance using the whole dataset. Alternatively, we used a robust version of standardization, which replaces mean and variance with median and median absolute deviation, correspondingly, and is applied on a per-plate basis (“Robust MAD”). Feature selection included variance thresholding to remove features with minimal variation across the dataset, removing highly correlated features, removing features with missing values or outliers, and removing blocklisted features—all using pycytominer32 default parameters.

> Way, G. P. et al. Predicting cell health phenotypes using image-based morphology profiling. Mol. Biol. Cell 32, 995–1005 (2021)

### How to run

This experment contains 4 Jupyter notebooks that:

0. [Download and preprocess profiles](./0_download_preprocess_data.ipynb) by performing two types of normalization (standardize and Robust MAD) and feature selection using [pycytominer](https://github.com/cytomining/pycytominer).
1. [Assess plate and well position effect](./1_plate_well_position_effect.ipynb) by calculating mAP using different metadata-based profile groupings across two types of normalization.
2. [Assess CRISPR guide phenotypic activity](./2_phenotypic_activity.ipynb) by calculating mAP for replicate retrieval against negative controls.
3. [Asses staining channel, cell compartment, and feature type contributions](./3_channel_compartment_analyses.ipynb) to phenotypic activity assesement using mAP.
4. [Assess perturbation phenotypic consistency](./4_phenotypic_consistency.ipynb) by retrieving guides targetting the same genes.
