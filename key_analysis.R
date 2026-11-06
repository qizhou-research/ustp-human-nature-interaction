# =============================================================================
# key_analysis.R
# "Undergrounding wastewater treatment plants turns urban nuisances into amenities"
# Corresponding author: Shen Qu (squ@bit.edu.cn)
#
# Core analytic settings behind the paper's main results. Assumes a cleaned,
# quality-screened survey data frame `data_full` (N = 780) is in scope.
# Individual-level data are withheld for participant privacy (see README).
# =============================================================================
library(dplyr)
library(lavaan)
library(MatchIt)
library(cobalt)

# ---- 1. Analytic sample and propensity-score matching (Fig. 2b) -------------
# Residents within 2,000 m of their surveyed plant (N = 581). 1:2 nearest-
# neighbour matching on a logistic propensity score, caliper 0.15 SD, without
# replacement. USTP = 1 underground plant, 0 aboveground plant.
match_sample <- function(d, seed = 20260813) {
    set.seed(seed)
    m <- matchit(USTP ~ Age + factor(Gender) + tenure_years + green_walk_time_cat + dist_to_plant_m,
                 data = d, method = "nearest", distance = "glm",
                 ratio = 2, caliper = 0.15, replace = FALSE)
    print(bal.tab(m, stats = "mean.diffs", un = TRUE))   # all |SMD| < 0.1
    match.data(m)
}
matched <- match_sample(filter(data_full, dist_to_plant_m <= 2000))   # N = 509

# ---- 2. Structural equation model and effect decomposition (Figs. 3, 4a-b) --
covs <- "Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m"
sem_model <- gsub("COV", covs, '
  # Measurement model
  PE =~ PE1 + PE2 + PE3 + PE4          # personalized ecology
  BH =~ BH1 + BH2 + BH3 + BH4          # benefits for humans (well-being)
  MT =~ MT1 + MT2 + MT3 + MT4          # motivation
  OP =~ OP1 + OP2 + OP3 + OP4          # opportunity
  CP =~ CP1 + CP2 + CP3 + CP4          # capability
  PB =~ PB1 + PB2 + PB4                # perceived benefit
  KN =~ KN1 + KN2 + KN3 + KN4          # knowledge
  TR =~ TR1 + TR2 + TR3 + TR4          # trust
  Risk =~ PB + KN + TR                 # second-order risk perception
  IA_neg =~ r1*IA1 + r1*IA2            # sensory nuisance (push)
  IA_pos =~ r2*IA3 + r2*IA4            # facility attraction (pull)
  IA1 ~~ en1*IA1; IA2 ~~ en2*IA2       # labelled residual variances (held equal
  IA3 ~~ ep1*IA3; IA4 ~~ ep2*IA4       #   across groups in the multi-group models)

  # Structural model, covariates on every endogenous variable
  IA_neg ~ a1*USTP + COV
  IA_pos ~ a2*USTP + COV
  Risk ~ b1*IA_neg + b2*IA_pos + COV
  CP   ~ c1*IA_neg + c2*IA_pos + COV
  OP   ~ d1*IA_neg + d2*IA_pos + COV
  MT   ~ e1*OP + e2*CP + e3*Risk + f1*IA_neg + f2*IA_pos + COV
  PE   ~ g1*OP + g2*CP + g3*MT + g4*Risk + h1*IA_neg + h2*IA_pos + COV
  BH   ~ i1*PE + j1*IA_neg + j2*IA_pos + k1*USTP + COV

  # Decomposition of the USTP effect (delta-method CIs)
  to_PE_neg := h1 + b1*g4 + b1*e3*g3 + c1*g2 + c1*e2*g3 + d1*g1 + d1*e1*g3 + f1*g3
  to_PE_pos := h2 + b2*g4 + b2*e3*g3 + c2*g2 + c2*e2*g3 + d2*g1 + d2*e1*g3 + f2*g3
  PE_nuisance   := a1 * to_PE_neg
  PE_attraction := a2 * to_PE_pos
  PE_total      := PE_nuisance + PE_attraction
  BH_nuisance   := a1 * (j1 + to_PE_neg*i1)
  BH_attraction := a2 * (j2 + to_PE_pos*i1)
  BH_direct     := k1
  BH_total      := BH_nuisance + BH_attraction + BH_direct
')
fit_sem <- function(model, d, ...)
    sem(model, data = d, estimator = "MLR", std.lv = TRUE, check.lv.names = FALSE, ...)

fit <- fit_sem(sem_model, matched)
print(fitMeasures(fit, c("chisq.scaled", "df.scaled", "cfi.robust", "rmsea.robust", "srmr")))
decomp <- function(f) filter(parameterEstimates(f), op == ":=", !grepl("^to_", lhs))[,
                             c("lhs", "est", "se", "pvalue", "ci.lower", "ci.upper")]
print(decomp(fit))

# ---- 3. Distance-threshold sensitivity (Fig. 4c, Supplementary Table S10) ---
for (th in c(1000, 1500, 2500, 3000, Inf)) {
    d <- match_sample(filter(data_full, dist_to_plant_m <= th))
    cat("\nThreshold", th, "m, matched N =", nrow(d), "\n")
    print(decomp(fit_sem(sem_model, d)))
}

# ---- 4. Subgroup heterogeneity, full-model multi-group SEM (Fig. 5) ---------
# Every labelled structural path is freed across the two groups and compared
# against a fully constrained model (likelihood-ratio test, df = 23). The free
# model has local optima that depend on group order, so it is fitted with both
# orders and the higher-likelihood solution is kept. Group names are read back
# from lavaan, which orders groups as given in group.label.
paths <- c("a1","a2","b1","b2","c1","c2","d1","d2","e1","e2","e3",
           "f1","f2","g1","g2","g3","g4","h1","h2","i1","j1","j2","k1")
base_model <- paste(grep(":=", strsplit(sem_model, "\n")[[1]], value = TRUE, invert = TRUE),
                    collapse = "\n")
multigroup <- function(group, drop_cov) {
    m0 <- gsub(paste0("\\s*\\+\\s*", drop_cov, "\\b"), "", base_model, perl = TRUE)
    m1 <- m0
    for (p in paths) m1 <- gsub(paste0("\\b", p, "\\*"), sprintf("c(%s_g1, %s_g2)*", p, p), m1, perl = TRUE)
    levs <- levels(matched[[group]])
    f0 <- fit_sem(m0, matched, group = group, group.label = levs)
    fits <- lapply(list(levs, rev(levs)), function(l) fit_sem(m1, matched, group = group, group.label = l))
    f1 <- fits[[which.max(sapply(fits, function(f) as.numeric(logLik(f))))]]
    print(anova(f1, f0))
    out <- filter(standardizedSolution(f1, ci = TRUE), op == "~", grepl("_g[12]$", label))[,
                  c("group", "lhs", "rhs", "est.std", "ci.lower", "ci.upper", "pvalue")]
    out$group <- lavInspect(f1, "group.label")[out$group]
    out
}
matched <- matched %>% mutate(
    GenderGroup = factor(ifelse(Gender == 1, "Male", "Female"), levels = c("Male", "Female")),
    AgeGroup    = factor(ifelse(Age > 45, "Older", "Younger"), levels = c("Older", "Younger")),
    TenureGroup = factor(ifelse(tenure_years > median(tenure_years, na.rm = TRUE), "Long", "Short"),
                         levels = c("Long", "Short")))
het_gender <- multigroup("GenderGroup", "Gender")
het_age    <- multigroup("AgeGroup", "Age")
het_tenure <- multigroup("TenureGroup", "tenure_years")
