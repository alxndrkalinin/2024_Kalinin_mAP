## 2. Cell Health

This experiment demostrated versatility of the mAP framework through its application to different tasks on “Cell Health” dataset, evaluating the effects of selected preprocessing methods and experimental designs.

### Dataset description

We used our previously published “Cell Health” dataset of Cell Painting images of CRISPR-Cas9 knockout perturbations of 59 genes, targeted by 119 guides in three different cell lines (A549, ES2, and HCC44). Morphological profiles were previously extracted from images using CellProfiler47 and median-aggregated on the well level. We used a subset of 100 guides that had exactly six replicates (two replicates in three different plates) in each cell line. We performed two types of profile normalization followed by feature selection using Pycytominer. The first preprocessing method included data standardization by subtracting means from feature values and dividing them by variance using the whole dataset. Alternatively, we used a robust version of standardization, which replaces mean and variance with median and median absolute deviation, correspondingly, and is applied on a per-plate basis (“Robust MAD”). Feature selection included variance thresholding to remove features with minimal variation across the dataset, removing highly correlated features, removing features with missing values or outliers, and removing blocklisted features—all using Pycytominer default parameters.

> Way, G. P. et al. Predicting cell health phenotypes using image-based morphology profiling. Mol. Biol. Cell 32, 995–1005 (2021)

### How to run

To download data and execute all analyses, run:

```bash
bash run_all_cellhealth.sh
```

Or run individual steps:

0. Download data: [`bash 0_download_data.sh`](./0_download_data.sh)
1. Preprocess profiles: [`python 1_preprocess_data.py`](./1_preprocess_data.py)
2. Calculate plate and well position effects on phenotypic activity: [`python 2_plate_well_position_effect.py`](./2_plate_well_position_effect.py)
3. Calculate phenotypic activity per channel and compartment: [`python 3_channel_compartment_analysis.py`](./3_channel_compartment_analysis.py)
4. Calculate phenotypic activity for both types of controls: [`4_phenotypic_activity.py`](./4_phenotypic_activity.py)
5. Calculate phenotypic consistency: [`5_phenotypic_consistency.py`](./5_phenotypic_consistency.py)
