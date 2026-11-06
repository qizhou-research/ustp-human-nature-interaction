# ustp-human-nature-interaction

Key analysis code for "Undergrounding wastewater treatment plants turns urban nuisances into amenities" (corresponding author Shen Qu, squ@bit.edu.cn).

`key_analysis.R` (about 100 lines) contains the analytic settings behind the paper's main results.

| Section | What it does | Reported in |
|---|---|---|
| 1 | 2,000 m analytic sample and 1:2 propensity-score matching | Fig. 2b |
| 2 | Structural equation model and decomposition of the USTP effect into nuisance, attraction and direct components | Figs. 3 and 4a,b, Extended Data 2 |
| 3 | Re-matching and re-estimation at alternative distance thresholds | Fig. 4c, Supplementary Table S10 |
| 4 | Full-model multi-group SEM by gender, age and residential tenure | Fig. 5, Supplementary Fig. S4 |

The script assumes a cleaned survey data frame `data_full` (N = 780 valid responses) with item columns PE1-4, CP1-4, OP1-4, BH1-4, MT1-4, IA1-4, PB1, PB2, PB4, KN1-4, TR1-4, and the variables `USTP` (1 underground, 0 aboveground), `Age`, `Gender` (1 male, 0 female), `tenure_years`, `green_walk_time_cat` (walking time to the nearest green space, 4 categories) and `dist_to_plant_m` (GPS distance to the surveyed plant). Construct abbreviations follow the paper's Methods.

## Requirements

R 4.x with

```r
install.packages(c("dplyr", "lavaan", "MatchIt", "cobalt"))
```

## Data

Individual-level survey data include residential locations and are withheld to protect participant privacy. De-identified data are available from the corresponding author on reasonable request.

## License

MIT, see `LICENSE`.

## Citation

Please cite the paper. Archived on Zenodo, DOI: [to be added after release].
