## 1. Simulations

This experiment performed comparing multiple evaluation metrics on simulated profile data.

The code generated perturbation and control profiles such that each perturbation had 2 to 4 replicates, in experiments with 12, 24, or 36 control replicates (in all combinations). The simulated profiles varied in feature size, ranging from 100 to 5,000, representing a typical range for various kinds of profiles, particularly after feature reduction or selection. For perturbation replicates, a certain proportion of features (ranging from 1% to 64%) were sampled from a normal distribution with a shifted mean 𝒩(1,1), while the rest were drawn from 𝒩(0,1) (Figure 2). Following the previously described method, we calculated mAP scores and corresponding p-values at the perturbation level to assess phenotypic activity for each perturbation. We measured performance by calculating recall as the proportion of the 100 simulated perturbations for which each metric achieved statistical significance (p < 0.05). Figure S2-S4 explored using other distributions and/or parameters for perturbed features. Figure S5 compared mAP performance with different distance measures (cosine, Euclidean, correlation).

### How to run

0. Specify parameters in [`config.yaml`](./config.yaml) and run `python 0_run_simulations.py`.
1. Open [1_plot_simulations_results.ipynb](./1_plot_simulations_results.ipynb) to plot results.

### Figures S2-5

For reproducing results shown in supplementary figures, overrride the corresponding parameters when running the script, e.g.:

```
python 0_run_simulations.py perturb_distribs.gaussian=[2,2]
```

to reproduce results in Supplementary Figure S3.