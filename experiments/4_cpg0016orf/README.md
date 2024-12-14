## 4. cpg0016[orf] (aka cpg0016-jump[orf])

This experiment demostrated application of mAP to image-based morphological profiles of cells perturbed by gene overexpression.

### Dataset description

We used the [JUMP Consortium’s](https://broad.io/jump) “cpg0016-jump[orf]” dataset (abbreviated to “cpg0016[orf]” here), which contains Cell Painting images of U2OS cells treated with 15,136 overexpression reagents (ORFs) encompassing 12,602 unique genes. Morphological profiles were previously extracted from images using CellProfiler, mean-aggregated on the well level, and then corrected for plate layout effects by subtracting means from feature values per well location. Cell counts were regressed out from each feature with more than 100 unique values. After that, profiles were normalized per plate by subtracting medians from feature values and dividing them by median absolute deviation (“Robust MAD”). Feature selection was performed using pycytominer and profiles were corrected for batch effects by a combination of the sphering transformation and Harmony (an iterative algorithm based on expectation-maximization that alternates between finding clusters with high diversity of batches, and computing mixture-based corrections within such clusters).

> Chandrasekaran, S. N. et al. Morphological map of under- and over-expression of genes in human cells. bioRxiv 2024.12.02.624527 (2024)

### How to run

0. [Run snakemake pipeline](./Snakefile) that downloads profiles, preprocesses them, and calculates mAP for both phenotypic activity and consistency assesement:

```bash
snakemake -c1
```
where `1` is the number of cores to use.

1. [Plot phenotypic activity](./1_phenotypic_activity.ipynb)
2. [Plot phenotypic consistency](./2_phenotypic_consistency.ipynb)
