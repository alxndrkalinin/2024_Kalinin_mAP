## 6. Perturb-seq dataset

This experiment demostrated application of mAP to single-cell mRNA profiling data.

### Dataset description

We used the public Perturb-seqs mRNA profiling dataset of single cells treated with CRISPRi containing 10X single-cell gene expression reads, barcode identities, and activity readouts (Gene Expression Omnibus accession GSE132080). The experiment assessed how single-guide RNAs (sgRNAs) containing mismatches to their target sites attenuate expression levels of target genes45. Specifically, 25 genes involved in a diverse range of essential cell biological processes were targeted with 5–6 mismatched sgRNAs, covering the range from full to low activity, and 10 nontargeting controls. Each mismatched guide was characterized by its activity levels relative to the perfectly matched sgRNA targeting the same gene. The distributions of sgRNAs were largely unimodal, although broader than those with the perfectly matched sgRNA or the control sgRNA45. We performed single-cell profile normalization and feature selection using Seurat.

> Jost M, Santos DA, Saunders RA, Horlbeck MA et al. Titrating gene expression using libraries of systematically attenuated CRISPR guide RNAs. Nat Biotechnol 2020 Mar;38(3):355-364. PMID: 31932729

### How to run

First, create aa additional new conda environment for data preprocessing using R:

```
conda env create -f perturbseq_processing_environment.yml
```

To download data and execute all analyses, run:

```bash
bash run_all_perturbseq.sh
```

Or run individual steps:

0. Download data: [`bash 0_download_data.sh`](./0_download_data.sh)
1. Preprocess perturb-seq data using R: [`conda run -n perturbseq-processing Rscript 1_process_perturbseq.r`](./1_process_perturbseq.r)
2. Finalize preprocessing in Python: [`python 2_finalize_perturbseq.py`](./2_finalize_perturbseq.py)
3. Calculate mAP [`python 3_calculate_map_perturbseq.py`](./3_calculate_map_perturbseq.py)
4. Plot mAP results: [`python 4_plot_map_perturbseq.py`](./4_plot_map_perturbseq.py)
