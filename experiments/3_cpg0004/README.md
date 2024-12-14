## 3. cpg0004 (aka LINCS)

This experiment demostrated application of mAP to image-based morphological profiles of cells perturbed by small molecules.

### Dataset description

We used our previously published dataset “cpg0004-lincs” (abbreviated to cpg0004 here) that contains Cell Painting images of 1,327 small-molecule perturbations of A549 human cells. The wells on each plate were perturbed with 56 different compounds in six different doses. Every compound was replicated 4 times per dose, with each replicated on a different plate. In this study, only the highest dose point of 10 μM was used. Morphological profiles were previously extracted from images using CellProfiler. Profile normalization, feature selection, and batch correction were performed using pycytominer. First, profiles were normalized against DMSO controls by subtracting medians from DMSO feature values and dividing them by median absolute deviation (“Robust MAD”). Feature selection included variance thresholding to remove features with minimal variation across the dataset, removing highly correlated features, and removing features with missing values or outliers. Finally, profiles were corrected for batch effects by the sphering transformation (computes a whitening transformation matrix based on negative controls and applies this transformation to the entire dataset).

> Way, G. P. et al. Morphology and gene expression profiling provide complementary information for mapping cell state. Cell Syst 13, 911–923.e9 (2022)

### How to run

0. [Run snakemake pipeline](./Snakefile) that downloads profiles, preprocesses them, and calculates mAP for both phenotypic activity and consistency assesement:

```bash
snakemake -c1
```
where `1` is the number of cores to use.

1. [Plot phenotypic activity](./1_phenotypic_activity.ipynb)
2. [Plot phenotypic consistency](./2_phenotypic_consistency.ipynb)
