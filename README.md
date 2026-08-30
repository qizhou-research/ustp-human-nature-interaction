# Analysis code

Code accompanying "Undergrounding sewage treatment reshapes human-nature interaction through positive ecological amenities rather than passive nuisance mitigation." Corresponding author Shen Qu (squ@bit.edu.cn).

This repository contains `analysis.R`, which documents the analytic settings behind the manuscript's results. It shows the propensity-score-matching specification, every structural equation model specification (the main model, its mediation decomposition, subgroup heterogeneity by gender/age/tenure, all robustness and sensitivity variants reported in the manuscript and its Supplementary Information, and the measurement-model checks behind Supplementary Tables S4 through S8 and Extended Data 3), and the fit-index criteria used to evaluate them.

`analysis.R` is a reference for what was fitted and how, not an executable end-to-end pipeline. Raw survey data are individual-level and are withheld from this public repository to protect participant privacy, consistent with the manuscript's Data availability statement, so the script assumes a cleaned analytic data frame is already in scope rather than reconstructing it from the original survey exports (see the comment at the top of section 1 in `analysis.R` for the assumed data structure). Figure-generation code and one-off GIS/geocoding data preparation are also not included, since neither documents an analytic setting.

## Requirements

The model-fitting calls in `analysis.R` use these R packages.

```r
install.packages(c("dplyr", "tidyr", "tibble", "psych", "lavaan", "semTools", "MatchIt", "cobalt"))
```

## How to read this file

`analysis.R` is organized to mirror the manuscript and Supplementary Information, with a comment above each block naming the table or figure it supports.

1. Data preparation (assumed clean sample, main analytic threshold)
2. Propensity-score matching specification
3. Main structural equation model
4. Mediation decomposition (Fig. 3)
5. Subgroup heterogeneity (gender, age, tenure)
6. Full-model multi-group SEM by subgroup (Fig. 4, Supplementary Fig. S4)
7. Robustness and sensitivity checks (Supplementary Tables S10 through S13, S16)
8. Measurement-model equality constraints (Supplementary Table S16)
9. Odor-perceptibility and odor-source correlations (Discussion)
10. Descriptive statistics, reliability, and CFA behind Supplementary Tables S4 through S8 and Extended Data 3 (basic demographic and item-level comparison tables, Tables S2-S3, are omitted as routine descriptive computations)

Each `lavaan` model string is the actual specification used to produce the corresponding reported numbers. Fit indices are reported on the same MLR-robust scale used throughout the manuscript's Methods, "Model fit indices" (chi-square to degrees of freedom ratio, CFI, RMSEA, SRMR).

## Data

Individual-level survey data, including the raw responses and the manual plant-assignment corrections used during cleaning, are available from the corresponding author upon reasonable request, consistent with the manuscript's Data availability statement.

## License

Code is released under the MIT License. See `LICENSE`.

## Citation

If you use this code, please cite the paper. This repository is archived on Zenodo at the following DOI, which resolves to the exact v1.0.0 release.

▎ DOI: 10.5281/zenodo.22178919 (https://doi.org/10.5281/zenodo.22178919)

Repository: https://github.com/qizhou-research/ustp-human-nature-interaction
