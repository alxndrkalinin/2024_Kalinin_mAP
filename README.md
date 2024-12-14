# A versatile information retrieval framework for evaluating profile strength and similarity

This repository contains the source code to reproduce the results in the
paper: "[A versatile information retrieval framework for evaluating profile strength and similarity](https://doi.org/10.1101/2024.04.01.587631)".

## Getting started

### System requirements
copairs supports Python 3.9+ and should work with all modern operating systems (tested with MacOS 13.5, Ubuntu 18.04).

### Dependencies
This code depends on widely used Python packages:
* numpy
* scipy
* pandas
* jupyter
* seaborn
* networkx
* umap-learn
* scikit-learn

It also uses [pycytominer](https://github.com/alxndrkalinin/pycytominer/tree/fix-ops-custom-features) for profiling data preprocessing and [copairs](https://github.com/cytomining/copairs/tree/v0.4.2) for profile grouping and mAP calculations.

### Installation

We suggest using [Conda](https://docs.conda.io/projects/conda/en/stable/) for
environment management. The following commands create the environment from
scratch and install the required packages.

```bash
conda create -n map_eval "python>=3.9"
conda activate map_eval
pip install .
```

## Contents

Results are organized per dataset in the [experiments](./experiments/) subdirectory.

Each experiment directory  includes brief description of the dataset and scripts and/or Jupyter notebooks to download and preprocess data, calculate metrics, and generate figures for the paper.

### Experiments

1. [Simulations](./experiments/1_simulations/) (Figures 2, S2-5)
2. [CellHealth data](./experiments/2_cellhealth/) (Figures 3, S6)
3. [cpg0004 data](./experiments/3_cpg0004/) (Figures S7A, S7C)
4. [cpg0016orf data](./experiments/4_cpg0016orf/) (Figures S7B, S7D)
5. [nELISA data](./experiments/5_nelisa/) (Figures 4A-B)
6. [PerturbSeq data](./experiments/6_perturbseq/) (Figures 4C-D, S5A-B, S8)
7. [Mitocheck data](./experiments/7_mitocheck/) (Figure 5C-D, S9-10)

## Citation

```bibtex
@article {Kalinin2024.04.01.587631,
        author = {Kalinin, Alexandr A. and Arevalo, John and Vulliard, Loan and Serrano, Erik and Tsang, Hillary and Bornholdt, Michael and Rajwa, Bartek and Carpenter, Anne E. and Way, Gregory P. and Singh, Shantanu},
        title = {A versatile information retrieval framework for evaluating profile strength and similarity},
        elocation-id = {2024.04.01.587631},
        year = {2024},
        doi = {10.1101/2024.04.01.587631},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2024/04/02/2024.04.01.587631},
        journal = {bioRxiv}
}
```