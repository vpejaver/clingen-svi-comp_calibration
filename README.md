## Code to calibrate tools for clinical interpretation and generate summarized results

This repository contains only the **code** relevant to the paper. Due to large file sizes, data and intermediate result files are hosted [here](https://mountsinai.box.com/s/x9nlvdxaqgznfy6sn7fo3je4qz99huc8). Please refer to the README at this link for information on the data files.

### Repository structure

The repository is organized as follows:
```bash
├── LICENSE
├── README.md
├── local_posterior_probability
│   ├── get_all_thresholds.m
│   ├── get_both_bootstrapped_posteriors.m
│   ├── get_both_local_posteriors.m
│   ├── get_discounted_thresholds.m
│   ├── main.m
│   ├── plot_both_posteriors.m
│   └── print_thresholds.m
├── plotting
│   ├── plot_both_posteriors_pub.m
│   ├── plot_heatmap_gnomad_set.m
│   ├── plot_heatmap_lr_testset.m
│   └── plot_posterior_wrapper.m
└── results_postprocessing
    ├── assess_default_thresholds.m
    ├── calculate_coverage.m
    └── make_thr_table.m
```

#### 1. `local_posterior_probability`

This directory contains the actual implementation of the algorithm to calculate local posterior probabilities (as described in Figure 2 of the paper). The script `main.m` serves as the wrapper that calls all the other functions in this directory.


#### 2. `results_postprocessing`

This directory contains scripts to post-process outputs from `local_posterior_probability` and/or generate additional statistics and tables for the results.
* `make_thr_table.m` : script to generate and systematically print out the score thresholds in Table 2 (and Supplemental Table S1). Note that the format is not exactly as in the paper but it should be easy to update manually to align with the format in the paper.
* `assess_default_thresholds.m` : script to generate Table 3.
* `calculate_coverage.m` : script to generate Supplemental Table S2.


#### 3. `plotting`

This directory contains the code used to make the plots in the paper.
* `plot_posterior_wrapper.m` : wrapper script to plot Figure 3. This script calls the function `plot_both_posteriors_pub.m`, which generates each individual local posterior probability plot, i.e., the function is called 26 times for each of the 26 subplots inside Figure 3.
* `plot_both_posterior_pub.m` : function to plot a single publication-quality local posterior probability plot. Note that this more or less does the same thing that `plot_both_posteriors.m` in `local_posterior_probability` does but the resulting plot matches the look and feel of the ones in the paper. It is recommended that this function be used to visualize finalized results.
* `plot_heatmap_lr_testset.m` : script to plot the heatmap summarizing interval-based likelihood ratios on the validation set (Figure 4A).
* `plot_heatmap_gnomad_set.m` : script to plot the heatmap summarizing the fraction of gnomAD variants falling within each score interval (Figure 4B).
