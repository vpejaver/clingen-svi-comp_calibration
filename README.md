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

