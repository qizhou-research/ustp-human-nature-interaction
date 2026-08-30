# =============================================================================
# analysis.R
#
# "Undergrounding sewage treatment reshapes human-nature interaction through
#  positive ecological amenities rather than passive nuisance mitigation"
# Corresponding author: Shen Qu (squ@bit.edu.cn)
#
# This file documents the analytic settings behind the manuscript's results:
# the propensity-score-matching specification, every structural equation
# model specification (main model, mediation decomposition, subgroup
# heterogeneity, all robustness and sensitivity variants reported in the
# manuscript and its Supplementary Information, and the measurement-model
# checks behind Supplementary Tables S4-S8, and Extended Data 3), and the
# fit-index criteria used
# to evaluate them. It is a reference for what was fitted and how, not an
# executable end-to-end pipeline: raw-file ingestion and data cleaning are
# summarized rather than reproduced (see section 1), and it assumes a
# cleaned analytic data frame is already in scope.
#
# Deliberately NOT included (see README.md): all ggplot2 figure-generation
# code, one-off GIS/geocoding/land-use data preparation, and the igraph-based
# mediation_paths() cross-check (a redundant arithmetic sanity check on the
# lavaan ":=" decomposition below, not itself a source of any reported number).
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(psych)
library(lavaan)
library(semTools)
library(MatchIt)   # propensity-score matching
library(cobalt)    # matching-balance diagnostics (bal.tab)

# ---- Shared helper for the robustness/sensitivity model variants ----
# Every A1/A2/A4/A5 variant below follows the same pattern: fit a lavaan model,
# print its decomposed ":=" effects (if any), and report its MLR-robust fit
# indices (chisq/df, CFI, RMSEA, SRMR) on the same scale used throughout the
# paper (Methods, "Model fit indices.").
fit_and_report <- function(model_string, data, label,
                           decomp_labels = c("ind_IA_neg_path", "ind_IA_pos_path",
                                             "direct_USTP_BH", "total_USTP_BH"),
                           ...) {
    fit <- sem(model_string, data = data, estimator = "MLR", std.lv = TRUE,
              check.lv.names = FALSE, ...)
    cat(sprintf("\n==== %s ====\nconverged: %s\n", label, lavInspect(fit, "converged")))
    decomp <- NULL
    if (length(decomp_labels)) {
        decomp <- parameterEstimates(fit, standardized = TRUE) %>%
            dplyr::filter(op == ":=", lhs %in% decomp_labels) %>%
            dplyr::select(lhs, est, se, pvalue, ci.lower, ci.upper)
        print(decomp)
    }
    fitm <- fitMeasures(fit, c("chisq.scaled", "df.scaled", "cfi.robust", "rmsea.robust", "srmr"))
    cat("Model fit (chisq/df, CFI, RMSEA, SRMR):", round(fitm, 3), "\n")
    list(fit = fit, decomp = decomp, fitm = fitm)
}


# =============================================================================
# 1. Data preparation
# =============================================================================
# Assumes a cleaned, quality-screened survey sample is already in scope as
# `data_full` (N = 780; see Methods, "Data quality"), with each construct's
# item columns (PE1-4, CP1-4, OP1-4, BH1-4, MT1-4, IA1-4, PB1-2/PB4 [PB3
# dropped during cleaning], KN1-4, TR1-4), USTP/Type, Age, Gender (1 = male,
# 0 = female), tenure_years (residential tenure, years), and the GPS-derived
# distance to the surveyed plant (dist_to_plant_m, metres, not yet thresholded)
# already computed. Raw-file ingestion (individual-level survey responses and
# the manual community-to-plant correction file are withheld from this public
# repository to protect participant privacy; see README.md), quality
# screening, and the plant-location join that produce this sample from the
# original survey exports are summarized here rather than reproduced.

# Main analytic sample: residents within 2,000 m of their surveyed plant
# (sensitivity-tested against the full distance-cutoff sweep in section 7c
# below and Supplementary Table S10; the unrestricted sample loses
# significance on a1).
data <- data_full %>%
    dplyr::filter(dist_to_plant_m <= 2000)


# =============================================================================
# 2. Propensity-score matching
# =============================================================================
# 1:2 nearest-neighbour matching on a logistic propensity score, caliper = 0.15
# SD of the propensity score. dist_to_plant_m (GPS-derived, objective) is used
# instead of the survey's subjective "walking time to plant" item to avoid
# common-method bias with the psychological outcomes; Gender is included
# because pre-match group balance on Gender was borderline (P = .054).
m.out <- matchit(USTP ~ Age + green_walk_time_cat + dist_to_plant_m + tenure_years + factor(Gender),
                 data = data, method = "nearest", distance = "glm",
                 ratio = 2, caliper = 0.15, replace = FALSE)
matched_data <- match.data(m.out)

cat("\n---- PSM balance (standardized mean differences) ----\n")
print(bal.tab(m.out, stats = "mean.diffs", un = TRUE))


# =============================================================================
# 3. Main structural equation model
# =============================================================================
data.level1 <- matched_data

my_model <- '
      # 1. Measurement model
      PE   =~ PE1 + PE2 + PE3 + PE4
      BH   =~ BH1 + BH2 + BH3 + BH4
      MT   =~ MT1 + MT2 + MT3 + MT4
      OP   =~ OP1 + OP2 + OP3 + OP4
      CP   =~ CP1 + CP2 + CP3 + CP4

      PB   =~ PB1 + PB2 + PB4
      KN   =~ KN1 + KN2 + KN3 + KN4
      TR   =~ TR1 + TR2 + TR3 + TR4
      IA_neg   =~ r1*IA1 + r1*IA2
      IA_pos   =~ r2*IA3 + r2*IA4
      IA1 ~~ en1*IA1; IA2 ~~ en2*IA2
      IA3 ~~ ep1*IA3; IA4 ~~ ep2*IA4

      # 2. Second-order risk factor (loadings freely estimated, not equality-constrained)
      Risk =~ PB + KN + TR

      # 3. Structural model. Age/Gender/tenure_years/green_walk_time_cat/dist_to_plant_m are carried
      #    as covariates on every endogenous variable.
      IA_neg ~ a1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
      IA_pos ~ a2*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

      Risk   ~ b1*IA_neg + b2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
      CP     ~ c1*IA_neg + c2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
      OP     ~ d1*IA_neg + d2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

      MT     ~ e1*OP + e2*CP + e3*Risk + f1*IA_neg + f2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
      PE     ~ g1*OP + g2*CP + g3*MT + g4*Risk + h1*IA_neg + h2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
      BH     ~ i1*PE + j1*IA_neg + j2*IA_pos + k1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

      # 4. USTP -> BH mediation decomposition (analytic, delta-method CIs)
      raw_IA_neg_to_BH := j1 + h1*i1 + b1*g4*i1 + b1*e3*g3*i1 + c1*g2*i1 + c1*e2*g3*i1 + d1*g1*i1 + d1*e1*g3*i1 + f1*g3*i1
      raw_IA_pos_to_BH := j2 + h2*i1 + b2*g4*i1 + b2*e3*g3*i1 + c2*g2*i1 + c2*e2*g3*i1 + d2*g1*i1 + d2*e1*g3*i1 + f2*g3*i1
      ind_IA_neg_path := a1 * raw_IA_neg_to_BH
      ind_IA_pos_path := a2 * raw_IA_pos_to_BH
      direct_USTP_BH := k1
      total_USTP_BH := ind_IA_neg_path + ind_IA_pos_path + k1

      # USTP -> PE decomposition (same logic, one step short of BH)
      raw_IA_neg_to_PE := h1 + b1*g4 + b1*e3*g3 + c1*g2 + c1*e2*g3 + d1*g1 + d1*e1*g3 + f1*g3
      raw_IA_pos_to_PE := h2 + b2*g4 + b2*e3*g3 + c2*g2 + c2*e2*g3 + d2*g1 + d2*e1*g3 + f2*g3
      ind_IA_neg_to_PE := a1 * raw_IA_neg_to_PE
      ind_IA_pos_to_PE := a2 * raw_IA_pos_to_PE
      total_USTP_PE := ind_IA_neg_to_PE + ind_IA_pos_to_PE
    '

sem.fit <- sem(my_model, data = data.level1, estimator = "MLR",
              std.lv = TRUE, check.lv.names = FALSE)

fit_summary <- fitMeasures(sem.fit, c("chisq.scaled","df.scaled","pvalue.scaled",
                                      "cfi.robust","rmsea.robust","srmr"))
cat("\n---- Main model fit indices (Supplementary Table S8) ----\n")
print(round(fit_summary, 4))

cat("\n---- Main model full summary ----\n")
print(summary(sem.fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE, ci = TRUE))


# =============================================================================
# 4. Mediation decomposition (Fig. 3)
# =============================================================================
extract_mediation <- function(fit) {
    pe <- parameterEstimates(fit, standardized = TRUE) %>% dplyr::filter(op == ":=")
    out <- data.frame(estimate = pe$est, se = pe$se, pvalue = pe$pvalue,
                      ci_low = pe$ci.lower, ci_high = pe$ci.upper)
    rownames(out) <- pe$lhs
    out
}

med <- extract_mediation(sem.fit)
cat("\n---- USTP -> BH / USTP -> PE decomposition (Fig. 3) ----\n")
print(round(med, 4))

med_df <- data.frame(
    component = c("Total effect", "Direct effect", "Facility attraction pathway", "Sensory nuisance pathway"),
    estimate  = c(med["total_USTP_BH","estimate"], med["direct_USTP_BH","estimate"],
                  med["ind_IA_pos_path","estimate"], med["ind_IA_neg_path","estimate"]),
    ci_low    = c(med["total_USTP_BH","ci_low"], med["direct_USTP_BH","ci_low"],
                  med["ind_IA_pos_path","ci_low"], med["ind_IA_neg_path","ci_low"]),
    ci_high   = c(med["total_USTP_BH","ci_high"], med["direct_USTP_BH","ci_high"],
                  med["ind_IA_pos_path","ci_high"], med["ind_IA_neg_path","ci_high"])
)
cat("\nFig. 3b (effect on benefits for humans):\n"); print(med_df)

pe_df <- data.frame(
    component = c("Total effect", "Facility attraction pathway", "Sensory nuisance pathway"),
    estimate  = c(med["total_USTP_PE","estimate"], med["ind_IA_pos_to_PE","estimate"], med["ind_IA_neg_to_PE","estimate"]),
    ci_low    = c(med["total_USTP_PE","ci_low"], med["ind_IA_pos_to_PE","ci_low"], med["ind_IA_neg_to_PE","ci_low"]),
    ci_high   = c(med["total_USTP_PE","ci_high"], med["ind_IA_pos_to_PE","ci_high"], med["ind_IA_neg_to_PE","ci_high"])
)
cat("\nFig. 3a (effect on personalized ecology):\n"); print(pe_df)


# =============================================================================
# 5. Subgroup heterogeneity (gender, age, residential tenure)
# =============================================================================

# ---- 5a. Age: lightweight multi-group SEM ----
data.level1$AgeGroup <- factor(ifelse(matched_data$Age > 45, "Older", "Younger"), levels = c("Older","Younger"))
mg_model_free <- '
        PE =~ PE1+PE2+PE3+PE4
        BH =~ BH1+BH2+BH3+BH4
        MT =~ MT1+MT2+MT3+MT4
        PE ~ c(mt2pe_o, mt2pe_y)*MT
        BH ~ c(pe2bh_o, pe2bh_y)*PE
    '
mg_model_constrained <- '
        PE =~ PE1+PE2+PE3+PE4
        BH =~ BH1+BH2+BH3+BH4
        MT =~ MT1+MT2+MT3+MT4
        PE ~ MT
        BH ~ PE
    '
fit_age_free <- sem(mg_model_free, data = data.level1, group = "AgeGroup", std.lv = TRUE, estimator = "MLR")
fit_age_constrained <- sem(mg_model_constrained, data = data.level1, group = "AgeGroup",
                           group.equal = c("regressions"), std.lv = TRUE, estimator = "MLR")
cat("\n--- Multi-group SEM (Age >45 vs <=45): MT->PE, PE->BH free vs. constrained equal ---\n")
print(anova(fit_age_free, fit_age_constrained))
par_age <- parameterEstimates(fit_age_free) %>% dplyr::filter(label %in% c("mt2pe_o","mt2pe_y","pe2bh_o","pe2bh_y"))
print(par_age[, c("label","est","pvalue")])

# ---- 5b. Gender: multi-group SEM on a1 (USTP->IA_neg) / a2 (USTP->IA_pos) ----
data.level1$GenderGroup <- factor(ifelse(matched_data$Gender == 1, "Male", "Female"), levels = c("Male","Female"))
cat("\nGenderGroup N:\n"); print(table(data.level1$GenderGroup))
mg_gender_free <- '
        IA_neg =~ IA1 + IA2
        IA_pos =~ IA3 + IA4
        IA_neg ~ c(a1_male, a1_female)*USTP
        IA_pos ~ c(a2_male, a2_female)*USTP
    '
mg_gender_constrained <- '
        IA_neg =~ IA1 + IA2
        IA_pos =~ IA3 + IA4
        IA_neg ~ USTP
        IA_pos ~ USTP
    '
fit_gender_free <- sem(mg_gender_free, data = data.level1, group = "GenderGroup", std.lv = TRUE, estimator = "MLR")
fit_gender_constrained <- sem(mg_gender_constrained, data = data.level1, group = "GenderGroup",
                              group.equal = c("regressions"), std.lv = TRUE, estimator = "MLR")
cat("\n--- Multi-group SEM (Male vs Female): a1/a2 free vs. constrained equal ---\n")
print(anova(fit_gender_free, fit_gender_constrained))
par_gender <- parameterEstimates(fit_gender_free, standardized = TRUE, ci = TRUE) %>%
    dplyr::filter(label %in% c("a1_male","a1_female","a2_male","a2_female"))
print(par_gender[, c("lhs","rhs","label","std.all","ci.lower","ci.upper","pvalue")])

# ---- 5c. Residential tenure: same lightweight multi-group SEM as age ----
tenure_median <- median(matched_data$tenure_years, na.rm = TRUE)
cat("\nResidential tenure median (years):", tenure_median, "\n")
data.level1$TenureGroup <- factor(ifelse(matched_data$tenure_years > tenure_median, "Long", "Short"), levels = c("Long","Short"))
mg_model_long_free <- '
        PE =~ PE1+PE2+PE3+PE4
        BH =~ BH1+BH2+BH3+BH4
        MT =~ MT1+MT2+MT3+MT4
        PE ~ c(mt2pe_l, mt2pe_s)*MT
        BH ~ c(pe2bh_l, pe2bh_s)*PE
    '
mg_model_long_constrained <- '
        PE =~ PE1+PE2+PE3+PE4
        BH =~ BH1+BH2+BH3+BH4
        MT =~ MT1+MT2+MT3+MT4
        PE ~ MT
        BH ~ PE
    '
fit_long_lite_constrained <- sem(mg_model_long_constrained, data = data.level1, group = "TenureGroup",
                                 group.equal = c("regressions"), std.lv = TRUE, estimator = "MLR")
# The unconstrained model can converge to a markedly worse local optimum from
# lavaan's default start; starting from the constrained solution (cross-checked
# against an independent start = "simple" run, same result) finds the true optimum.
fit_long_lite_free <- sem(mg_model_long_free, data = data.level1, group = "TenureGroup",
                          std.lv = TRUE, estimator = "MLR", start = fit_long_lite_constrained)
cat("\n--- Multi-group SEM (Long vs Short tenure): MT->PE, PE->BH free vs. constrained equal ---\n")
print(anova(fit_long_lite_free, fit_long_lite_constrained))
par_long_lite <- parameterEstimates(fit_long_lite_free) %>%
    dplyr::filter(label %in% c("mt2pe_l","mt2pe_s","pe2bh_l","pe2bh_s"))
print(par_long_lite[, c("label","est","pvalue")])

# =============================================================================
# 6. Full-model (23-path) multi-group SEM by gender / age / tenure
#    (Fig. 4 and Supplementary Fig. S4; source of the reported subgroup path
#    coefficients. Supplementary Fig. S3, a per-subgroup mediation
#    decomposition, is a distinct, smaller analysis not computed here.)
# =============================================================================
strip_defs <- function(text) {
    lines <- strsplit(text, "\n")[[1]]
    paste(lines[!grepl(":=", lines)], collapse = "\n")
}
protected_labels <- c("a1","a2","b1","b2","c1","c2","d1","d2","e1","e2","e3",
                      "f1","f2","g1","g2","g3","g4","h1","h2","i1","j1","j2","k1")
make_group_free_structural <- function(text, s1, s2) {
    for (lab in protected_labels) {
        pattern <- paste0("\\b", lab, "\\*")
        repl <- paste0("c(", lab, "_", s1, ", ", lab, "_", s2, ")*")
        text <- gsub(pattern, repl, text, perl = TRUE)
    }
    text
}
my_model_base <- strip_defs(my_model)

run_full_multigroup_ci <- function(group_var, drop_covariate, s1, s2, label) {
    base_text <- gsub(paste0("\\s*\\+\\s*", drop_covariate, "\\b"), "", my_model_base, perl = TRUE)
    free_text <- make_group_free_structural(base_text, s1, s2)
    cat(sprintf("\n\n==== Full-model multi-group SEM (with CI): %s ====\n", label))
    fit_constrained <- sem(base_text, data = data.level1, group = group_var,
                           std.lv = TRUE, estimator = "MLR", check.lv.names = FALSE)
    fit_free <- sem(free_text, data = data.level1, group = group_var,
                    std.lv = TRUE, estimator = "MLR", check.lv.names = FALSE)
    cat("converged:", lavInspect(fit_free, "converged"),
        " admissible:", lavInspect(fit_free, "post.check"), "\n")
    omni <- anova(fit_free, fit_constrained)
    par_full <- standardizedSolution(fit_free, type = "std.all", ci = TRUE) %>%
        dplyr::filter(op == "~", grepl(paste0("_(", s1, "|", s2, ")$"), label)) %>%
        dplyr::select(group, lhs, rhs, label, est = est.std, ci.lower, ci.upper, pvalue)
    print(par_full)
    list(fit_free = fit_free, fit_constrained = fit_constrained, par = par_full,
         omni_chisq = omni[2, "Chisq diff"], omni_p = omni[2, "Pr(>Chisq)"])
}

par_pooled <- standardizedSolution(sem.fit, type = "std.all", ci = TRUE) %>%
    dplyr::filter(op == "~", label %in% protected_labels) %>%
    dplyr::select(label_base = label, est = est.std, ci.lower, ci.upper, pvalue)
cat("\n---- Pooled (main model) path estimates, for comparison ----\n"); print(par_pooled)

res_full_gender_ci <- run_full_multigroup_ci("GenderGroup", "Gender", "male", "female", "Male vs Female")
res_full_age_ci    <- run_full_multigroup_ci("AgeGroup", "Age", "older", "younger", "Older(>45) vs Younger(<=45)")
res_full_tenure      <- run_full_multigroup_ci("TenureGroup", "tenure_years", "long", "short", "Long-tenure vs Short-tenure")


# =============================================================================
# 7. Robustness and sensitivity checks (Supplementary Tables S10-S13, S16)
# =============================================================================

# ---- 7a. Risk operationalization (Supplementary Table S13, rows 1-4) ----
build_model_risk_single <- function(risk_indicators) {
    paste0('
          PE   =~ PE1 + PE2 + PE3 + PE4
          BH   =~ BH1 + BH2 + BH3 + BH4
          MT   =~ MT1 + MT2 + MT3 + MT4
          OP   =~ OP1 + OP2 + OP3 + OP4
          CP   =~ CP1 + CP2 + CP3 + CP4
          Risk =~ ', risk_indicators, '
          IA_neg   =~ r1*IA1 + r1*IA2
          IA_pos   =~ r2*IA3 + r2*IA4
          IA1 ~~ en1*IA1; IA2 ~~ en2*IA2
          IA3 ~~ ep1*IA3; IA4 ~~ ep2*IA4

          IA_neg ~ a1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          IA_pos ~ a2*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          Risk   ~ b1*IA_neg + b2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          CP     ~ c1*IA_neg + c2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          OP     ~ d1*IA_neg + d2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          MT     ~ e1*OP + e2*CP + e3*Risk + f1*IA_neg + f2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          PE     ~ g1*OP + g2*CP + g3*MT + g4*Risk + h1*IA_neg + h2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          BH     ~ i1*PE + j1*IA_neg + j2*IA_pos + k1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          raw_IA_neg_to_BH := j1 + h1*i1 + b1*g4*i1 + b1*e3*g3*i1 + c1*g2*i1 + c1*e2*g3*i1 + d1*g1*i1 + d1*e1*g3*i1 + f1*g3*i1
          raw_IA_pos_to_BH := j2 + h2*i1 + b2*g4*i1 + b2*e3*g3*i1 + c2*g2*i1 + c2*e2*g3*i1 + d2*g1*i1 + d2*e1*g3*i1 + f2*g3*i1
          ind_IA_neg_path := a1 * raw_IA_neg_to_BH
          ind_IA_pos_path := a2 * raw_IA_pos_to_BH
          direct_USTP_BH := k1
          total_USTP_BH := ind_IA_neg_path + ind_IA_pos_path + k1
        ')
}
fit_risk_trust <- fit_and_report(build_model_risk_single("TR1 + TR2 + TR3 + TR4"), data.level1, "Risk = Trust-only")$fit
fit_risk_pb    <- fit_and_report(build_model_risk_single("PB1 + PB2 + PB4"), data.level1, "Risk = Perceived-benefit-only")$fit
fit_risk_kn    <- fit_and_report(build_model_risk_single("KN1 + KN2 + KN3 + KN4"), data.level1, "Risk = Knowledge-only")$fit

model_risk_firstorder <- '
          PE   =~ PE1 + PE2 + PE3 + PE4
          BH   =~ BH1 + BH2 + BH3 + BH4
          MT   =~ MT1 + MT2 + MT3 + MT4
          OP   =~ OP1 + OP2 + OP3 + OP4
          CP   =~ CP1 + CP2 + CP3 + CP4
          PB   =~ PB1 + PB2 + PB4
          KN   =~ KN1 + KN2 + KN3 + KN4
          TR   =~ TR1 + TR2 + TR3 + TR4
          IA_neg   =~ r1*IA1 + r1*IA2
          IA_pos   =~ r2*IA3 + r2*IA4
          IA1 ~~ en1*IA1; IA2 ~~ en2*IA2
          IA3 ~~ ep1*IA3; IA4 ~~ ep2*IA4

          IA_neg ~ a1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          IA_pos ~ a2*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          PB     ~ b1_PB*IA_neg + b2_PB*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          KN     ~ b1_KN*IA_neg + b2_KN*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          TR     ~ b1_TR*IA_neg + b2_TR*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          CP     ~ c1*IA_neg + c2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          OP     ~ d1*IA_neg + d2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          MT     ~ e1*OP + e2*CP + e3_PB*PB + e3_KN*KN + e3_TR*TR + f1*IA_neg + f2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          PE     ~ g1*OP + g2*CP + g3*MT + g4_PB*PB + g4_KN*KN + g4_TR*TR + h1*IA_neg + h2*IA_pos + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          BH     ~ i1*PE + j1*IA_neg + j2*IA_pos + k1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          raw_IA_neg_to_PE := h1 + (b1_PB*g4_PB + b1_KN*g4_KN + b1_TR*g4_TR) +
                              (b1_PB*e3_PB*g3 + b1_KN*e3_KN*g3 + b1_TR*e3_TR*g3) +
                              c1*g2 + c1*e2*g3 + d1*g1 + d1*e1*g3 + f1*g3
          raw_IA_pos_to_PE := h2 + (b2_PB*g4_PB + b2_KN*g4_KN + b2_TR*g4_TR) +
                              (b2_PB*e3_PB*g3 + b2_KN*e3_KN*g3 + b2_TR*e3_TR*g3) +
                              c2*g2 + c2*e2*g3 + d2*g1 + d2*e1*g3 + f2*g3
          raw_IA_neg_to_BH := j1 + raw_IA_neg_to_PE*i1
          raw_IA_pos_to_BH := j2 + raw_IA_pos_to_PE*i1
          ind_IA_neg_path := a1 * raw_IA_neg_to_BH
          ind_IA_pos_path := a2 * raw_IA_pos_to_BH
          direct_USTP_BH := k1
          total_USTP_BH := ind_IA_neg_path + ind_IA_pos_path + direct_USTP_BH
        '
fit_risk_firstorder <- fit_and_report(model_risk_firstorder, data.level1,
                                      "Risk first-order (PB/KN/TR not pooled)")$fit

# ---- 7b. IA_pos operationalization (Supplementary Table S12) ----
model_iapos_split <- '
          PE   =~ PE1 + PE2 + PE3 + PE4
          BH   =~ BH1 + BH2 + BH3 + BH4
          MT   =~ MT1 + MT2 + MT3 + MT4
          OP   =~ OP1 + OP2 + OP3 + OP4
          CP   =~ CP1 + CP2 + CP3 + CP4
          PB   =~ PB1 + PB2 + PB4
          KN   =~ KN1 + KN2 + KN3 + KN4
          TR   =~ TR1 + TR2 + TR3 + TR4
          Risk =~ PB + KN + TR
          IA_neg   =~ r1*IA1 + r1*IA2
          IA1 ~~ en1*IA1; IA2 ~~ en2*IA2

          IA_neg ~ a1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          IA3 ~ a2v*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          IA4 ~ a2e*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          Risk   ~ b1*IA_neg + b2v*IA3 + b2e*IA4 + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          CP     ~ c1*IA_neg + c2v*IA3 + c2e*IA4 + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          OP     ~ d1*IA_neg + d2v*IA3 + d2e*IA4 + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          MT     ~ e1*OP + e2*CP + e3*Risk + f1*IA_neg + f2v*IA3 + f2e*IA4 + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          PE     ~ g1*OP + g2*CP + g3*MT + g4*Risk + h1*IA_neg + h2v*IA3 + h2e*IA4 + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
          BH     ~ i1*PE + j1*IA_neg + j2v*IA3 + j2e*IA4 + k1*USTP + Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m

          raw_IA_neg_to_PE := h1 + b1*g4 + b1*e3*g3 + c1*g2 + c1*e2*g3 + d1*g1 + d1*e1*g3 + f1*g3
          raw_IA3_to_PE    := h2v + b2v*g4 + b2v*e3*g3 + c2v*g2 + c2v*e2*g3 + d2v*g1 + d2v*e1*g3 + f2v*g3
          raw_IA4_to_PE    := h2e + b2e*g4 + b2e*e3*g3 + c2e*g2 + c2e*e2*g3 + d2e*g1 + d2e*e1*g3 + f2e*g3
          raw_IA_neg_to_BH := j1  + raw_IA_neg_to_PE*i1
          raw_IA3_to_BH    := j2v + raw_IA3_to_PE*i1
          raw_IA4_to_BH    := j2e + raw_IA4_to_PE*i1
          ind_IA_neg_path := a1  * raw_IA_neg_to_BH
          ind_IA3_path    := a2v * raw_IA3_to_BH
          ind_IA4_path    := a2e * raw_IA4_to_BH
          direct_USTP_BH  := k1
          total_USTP_BH   := ind_IA_neg_path + ind_IA3_path + ind_IA4_path + direct_USTP_BH
          ind_pos_combined := ind_IA3_path + ind_IA4_path
        '
res_iapos_split <- fit_and_report(
    model_iapos_split, data.level1, "IA_pos split into IA3 (visitation) / IA4 (interaction quality)",
    decomp_labels = c("ind_IA_neg_path","ind_IA3_path","ind_IA4_path","direct_USTP_BH","total_USTP_BH","ind_pos_combined"))

# ---- 7c. Distance-cutoff sweep (Supplementary Table S10) ----
# Re-applies the same clean sample (data_full, section 1) at each of six
# distance thresholds, redoing propensity-score matching (same specification
# as section 2) and refitting the main model (same specification as section 3)
# at each one.
run_pipeline_at_threshold <- function(threshold, seed = 20260813) {
    d <- data_full %>%
        dplyr::filter(dist_to_plant_m <= threshold)

    set.seed(seed)
    m.out_t <- matchit(USTP ~ Age + green_walk_time_cat + dist_to_plant_m + tenure_years + factor(Gender),
                       data = d, method = "nearest", distance = "glm",
                       ratio = 2, caliper = 0.15, replace = FALSE)
    dl1_t <- match.data(m.out_t)

    res <- fit_and_report(my_model, dl1_t, sprintf("Distance threshold = %s", threshold),
                          decomp_labels = c("ind_IA_neg_path","ind_IA_pos_path","direct_USTP_BH","total_USTP_BH"))
    list(threshold = threshold, N = nrow(dl1_t), fitm = res$fitm)
}

threshold_grid <- c(1000, 1500, 2000, 2500, 3000, Inf)
threshold_label <- function(th) if (is.infinite(th)) "No cutoff" else paste0(th, "m")
cat("\n---- Distance-cutoff sweep (Supplementary Table S10) ----\n")
for (th in threshold_grid) {
    lbl <- threshold_label(th)
    if (th == 2000) {
        cat(sprintf("\n%s (main model, reused), N=%d\n", lbl, nrow(data.level1)))
        cat("total_USTP_BH:", round(med["total_USTP_BH","estimate"], 4), "\n")
        cat("ind_IA_pos_path:", round(med["ind_IA_pos_path","estimate"], 4), "\n")
        cat("fit chisq/df, CFI, RMSEA, SRMR:",
            round(fit_summary[c("chisq.scaled","df.scaled","cfi.robust","rmsea.robust","srmr")], 3), "\n")
    } else {
        run_pipeline_at_threshold(th)
    }
}

# ---- 7d. Covariate-set variants (Supplementary Table S13, rows 5-6) ----
build_model_covar <- function(covar_str) {
    gsub("%COVAR%", covar_str, '
          PE   =~ PE1 + PE2 + PE3 + PE4
          BH   =~ BH1 + BH2 + BH3 + BH4
          MT   =~ MT1 + MT2 + MT3 + MT4
          OP   =~ OP1 + OP2 + OP3 + OP4
          CP   =~ CP1 + CP2 + CP3 + CP4
          PB   =~ PB1 + PB2 + PB4
          KN   =~ KN1 + KN2 + KN3 + KN4
          TR   =~ TR1 + TR2 + TR3 + TR4
          IA_neg   =~ r1*IA1 + r1*IA2
          IA_pos   =~ r2*IA3 + r2*IA4
          IA1 ~~ en1*IA1; IA2 ~~ en2*IA2
          IA3 ~~ ep1*IA3; IA4 ~~ ep2*IA4
          Risk =~ PB + KN + TR

          IA_neg ~ a1*USTP %COVAR%
          IA_pos ~ a2*USTP %COVAR%
          Risk   ~ b1*IA_neg + b2*IA_pos %COVAR%
          CP     ~ c1*IA_neg + c2*IA_pos %COVAR%
          OP     ~ d1*IA_neg + d2*IA_pos %COVAR%
          MT     ~ e1*OP + e2*CP + e3*Risk + f1*IA_neg + f2*IA_pos %COVAR%
          PE     ~ g1*OP + g2*CP + g3*MT + g4*Risk + h1*IA_neg + h2*IA_pos %COVAR%
          BH     ~ i1*PE + j1*IA_neg + j2*IA_pos + k1*USTP %COVAR%

          raw_IA_neg_to_BH := j1 + h1*i1 + b1*g4*i1 + b1*e3*g3*i1 + c1*g2*i1 + c1*e2*g3*i1 + d1*g1*i1 + d1*e1*g3*i1 + f1*g3*i1
          raw_IA_pos_to_BH := j2 + h2*i1 + b2*g4*i1 + b2*e3*g3*i1 + c2*g2*i1 + c2*e2*g3*i1 + d2*g1*i1 + d2*e1*g3*i1 + f2*g3*i1
          ind_IA_neg_path := a1 * raw_IA_neg_to_BH
          ind_IA_pos_path := a2 * raw_IA_pos_to_BH
          direct_USTP_BH := k1
          total_USTP_BH := ind_IA_neg_path + ind_IA_pos_path + direct_USTP_BH
        ', fixed = TRUE)
}
covar_variants <- list(
    covar_plant_distance = "+ Age + Gender + tenure_years + green_walk_time_cat + plant_walk_time_cat",  # self-reported distance, Table S13 row 5
    covar_minimal        = "+ Age"  # age-only, Table S13 row 6
)
covar_fits <- list()
cat("\n---- Covariate-set variants (Supplementary Table S13, rows 5-6) ----\n")
for (nm in names(covar_variants)) {
    covar_fits[[nm]] <- fit_and_report(build_model_covar(covar_variants[[nm]]), data.level1,
                                       sprintf("%s (%s)", nm, covar_variants[[nm]]))$fit
}

# ---- 7e. Cross-specification reliability of the theoretical backbone (Fig. 5a) ----
backbone_paths <- tibble::tribble(
    ~label, ~from,    ~to,
    "b1",   "IA_neg", "Risk",  "b2",   "IA_pos", "Risk",
    "c1",   "IA_neg", "CP",    "c2",   "IA_pos", "CP",
    "d1",   "IA_neg", "OP",    "d2",   "IA_pos", "OP",
    "f1",   "IA_neg", "MT",    "f2",   "IA_pos", "MT",
    "e1",   "OP",     "MT",    "e2",   "CP",     "MT",    "e3", "Risk", "MT",
    "g1",   "OP",     "PE",    "g2",   "CP",     "PE",
    "g3",   "MT",     "PE",    "g4",   "Risk",   "PE"
)
backbone_fits <- list(
    main = sem.fit, risk_trust_only = fit_risk_trust, risk_pb_only = fit_risk_pb, risk_kn_only = fit_risk_kn,
    covar_plantdist = covar_fits[["covar_plant_distance"]], covar_minimal = covar_fits[["covar_minimal"]]
)
backbone_long <- dplyr::bind_rows(lapply(names(backbone_fits), function(nm) {
    pe <- parameterEstimates(backbone_fits[[nm]], standardized = TRUE) %>%
        dplyr::filter(op == "~", label %in% backbone_paths$label)
    pe$variant <- nm
    pe
}))
backbone_summary <- backbone_long %>%
    dplyr::inner_join(backbone_paths, by = "label") %>%
    dplyr::group_by(label, from, to) %>%
    dplyr::summarise(
        n_variant = dplyr::n(), n_sig = sum(pvalue < .05),
        n_pos = sum(std.all > 0), n_neg = sum(std.all < 0),
        sign_consistent = (n_pos == n_variant) | (n_neg == n_variant),
        mean_std = mean(std.all), .groups = "drop"
    ) %>%
    dplyr::mutate(robust = sign_consistent & (n_sig == n_variant))
cat("\n---- Backbone-path reliability across", length(backbone_fits), "specifications (Fig. 5a) ----\n")
print(as.data.frame(backbone_summary[, c("label","from","to","n_sig","n_variant","sign_consistent","robust","mean_std")]))


# =============================================================================
# 8. Measurement-model equality constraints search (Supplementary Table S16)
# =============================================================================
# Tests only whether the same-label loadings on IA_neg/IA_pos/Risk are
# required for identification (SI Text S2, check 5). Residual variances of
# IA_neg/IA_pos are always freely estimated (matching my_model) and are not
# part of this search. A separate backward-elimination search over the Age
# control paths was in the original analysis script but is NOT reported: SI
# Text S2 check 5 explicitly states "We did not additionally search over the
# covariate specification (for example, dropping individual Age control
# paths)," so that step is intentionally omitted here.
lrt_row <- function(fit_new, fit_base, label) {
    a <- anova(fit_base, fit_new)
    data.frame(model = label, d.chisq = round(a$`Chisq diff`[2], 3),
              d.df = a$`Df diff`[2], p.value = round(a$`Pr(>Chisq)`[2], 4))
}
structural_full <- '
        IA_neg ~ a1*USTP + age1*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        IA_pos ~ a2*USTP + age2*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        Risk   ~ b1*IA_neg + b2*IA_pos + age3*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        CP     ~ c1*IA_neg + c2*IA_pos + age4*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        OP     ~ d1*IA_neg + d2*IA_pos + age5*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        MT     ~ e1*OP + e2*CP + e3*Risk + f1*IA_neg + f2*IA_pos + age6*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        PE     ~ g1*OP + g2*CP + g3*MT + g4*Risk + h1*IA_neg + h2*IA_pos + age7*Age + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
        BH     ~ i1*PE + j1*IA_neg + j2*IA_pos + age8*Age + k1*USTP + Gender + tenure_years + green_walk_time_cat + dist_to_plant_m
    '
common_meas_b1 <- '
        PE   =~ PE1 + PE2 + PE3 + PE4
        BH   =~ BH1 + BH2 + BH3 + BH4
        MT   =~ MT1 + MT2 + MT3 + MT4
        OP   =~ OP1 + OP2 + OP3 + OP4
        CP   =~ CP1 + CP2 + CP3 + CP4
        PB   =~ PB1 + PB2 + PB4
        KN   =~ KN1 + KN2 + KN3 + KN4
        TR   =~ TR1 + TR2 + TR3 + TR4
    '
build_meas_block <- function(ia_neg = c("baseline","free_loading"),
                             ia_pos = c("baseline","free_loading"),
                             risk   = c("baseline","free_loading")) {
    ia_neg <- match.arg(ia_neg); ia_pos <- match.arg(ia_pos); risk <- match.arg(risk)
    ia_neg_ld <- switch(ia_neg, baseline = "r1*IA1 + r1*IA2", free_loading = "ln1*IA1 + ln2*IA2")
    ia_pos_ld <- switch(ia_pos, baseline = "r2*IA3 + r2*IA4", free_loading = "lp1*IA3 + lp2*IA4")
    risk_ld   <- switch(risk, baseline = "l1*PB + l1*KN + l1*TR", free_loading = "PB + KN + TR")
    paste0("IA_neg =~ ", ia_neg_ld, "\n", "IA_pos =~ ", ia_pos_ld, "\n",
           "IA1~~en1*IA1; IA2~~en2*IA2\n", "IA3~~ep1*IA3; IA4~~ep2*IA4\n", "Risk   =~ ", risk_ld, "\n")
}
build_full_model <- function(meas_block, drop_age = character(0)) {
    extra <- if (length(drop_age)) paste(paste0(drop_age, " == 0"), collapse = "\n") else ""
    paste(common_meas_b1, meas_block, structural_full, extra, sep = "\n")
}

fit0 <- sem(build_full_model(build_meas_block()), data = data.level1, std.lv = TRUE)

meas_variant_specs <- list(ia_neg_free_loading = list(ia_neg = "free_loading"),
                           ia_pos_free_loading = list(ia_pos = "free_loading"),
                           risk_free_loading   = list(risk   = "free_loading"))
fit_A <- lapply(meas_variant_specs, function(spec) sem(build_full_model(do.call(build_meas_block, spec)), data = data.level1, std.lv = TRUE))
lrt_A <- do.call(rbind, Map(function(f, nm) lrt_row(f, fit0, nm), fit_A, names(fit_A)))
cat("\n---- Step A: measurement-model equality constraints, LRT vs. baseline ----\n"); print(lrt_A)

pick_variant <- function(lrt_tab, candidates) {
    sub <- lrt_tab[lrt_tab$model %in% candidates & lrt_tab$p.value < 0.05, ]
    if (nrow(sub) == 0) return("baseline")
    sub$model[which.max(sub$d.chisq)]
}
best_spec <- list(
    ia_neg = sub("^ia_neg_", "", pick_variant(lrt_A, "ia_neg_free_loading")),
    ia_pos = sub("^ia_pos_", "", pick_variant(lrt_A, "ia_pos_free_loading")),
    risk   = if (pick_variant(lrt_A, "risk_free_loading") == "risk_free_loading") "free_loading" else "baseline"
)
cat("Step A selected measurement-model specification:\n"); print(best_spec)


# =============================================================================
# 9. Odor-perceptibility and odor-source correlations (Discussion)
# =============================================================================
# Questionnaire item cross-reference: odor_noticeability has 3 levels
# ("not noticeable"/"moderate"/"noticeable"); Source1-4 = odor source (multi-
# select, Source1=industrial activities, Source2=sewage-treatment facility
# itself, Source3=vehicle traffic, Source4=other), only answered by
# respondents who did not select "not noticeable" for odor_noticeability (skip logic).
d_c <- data.level1 %>%
    dplyr::mutate(IA_neg_composite = (IA1 + IA2) / 2,
                  odor_noticeability_num = dplyr::case_when(odor_noticeability == "not_noticeable" ~ 1, odor_noticeability == "moderate" ~ 2, odor_noticeability == "noticeable" ~ 3))

cat("\n---- Odor noticeability distribution ----\n"); print(table(d_c$odor_noticeability, useNA = "ifany"))
cat("Odor-source endorsement rate (denominator = respondents who answered this item):\n")
print(round(sapply(d_c[, c("Source1","Source2","Source3","Source4")], function(x) mean(x == 1, na.rm = TRUE)), 3))

cat("\n---- Odor noticeability vs. sensory-nuisance perception (IA_neg) ----\n")
ct1 <- cor.test(d_c$odor_noticeability_num, d_c$IA_neg_composite)
cat("odor_noticeability (1-3) vs IA_neg composite: r =", round(ct1$estimate, 3), ", p =", round(ct1$p.value, 4), "\n")

source_labels <- c(Source1 = "Industrial activities", Source2 = "Sewage treatment facilities",
                   Source3 = "Vehicle traffic", Source4 = "Other sources")
cat("\n---- Odor source vs. sensory-nuisance perception (IA_neg) ----\n")
for (src in c("Source1","Source2","Source3","Source4")) {
    ok <- stats::complete.cases(d_c[[src]], d_c$IA_neg_composite)
    ct <- stats::cor.test(d_c[[src]][ok], d_c$IA_neg_composite[ok])
    cat(src, "(", source_labels[src], ") vs IA_neg composite: r =", round(ct$estimate, 3),
        ", p =", round(ct$p.value, 4), ", n =", ct$parameter + 2, "\n")
}


# =============================================================================
# 10. Descriptive statistics, reliability, and CFA behind Supplementary
#     Tables S4-S8, and Extended Data 3 (full structural path coefficients)
# =============================================================================
# All computed on `data` (dist_to_plant_m <= 2000 m, N = 581, pre-PSM), matching
# Methods' "N = 581" description. PB3 was already dropped upstream in section 1
# (select(-PB3)), so item-level tables here cover 31 of the 32 questionnaire
# items.
item_cols <- c(paste0("PE",1:4), paste0("CP",1:4), paste0("OP",1:4), paste0("BH",1:4),
              paste0("MT",1:4), paste0("IA",1:4), "PB1","PB2","PB4", paste0("KN",1:4), paste0("TR",1:4))

# ---- Table S4: skewness/kurtosis + reliability ----
norm_tab <- psych::describe(data[, item_cols])[, c("mean","sd","skew","kurtosis")]
class(norm_tab) <- "data.frame"  # psych::describe() returns a "psych"-classed object whose
                                  # own print method chokes once a character column is added
norm_tab$item <- rownames(norm_tab); rownames(norm_tab) <- NULL
norm_tab[, c("mean","sd","skew","kurtosis")] <- round(norm_tab[, c("mean","sd","skew","kurtosis")], 3)
cat("\n---- Table S4: skewness/kurtosis ----\n"); print(norm_tab[, c("item","mean","sd","skew","kurtosis")])

cat("\n---- Mardia's multivariate normality test ----\n")
print(semTools::mardiaSkew(data[, item_cols])); print(semTools::mardiaKurtosis(data[, item_cols]))

constructs <- list(PE=paste0("PE",1:4), CP=paste0("CP",1:4), OP=paste0("OP",1:4), BH=paste0("BH",1:4),
                   MT=paste0("MT",1:4), IA_neg=c("IA1","IA2"), IA_pos=c("IA3","IA4"),
                   PB=c("PB1","PB2","PB4"), KN=paste0("KN",1:4), TR=paste0("TR",1:4))
alpha_tab <- do.call(rbind, lapply(names(constructs), function(nm) {
    a <- suppressWarnings(psych::alpha(data[, constructs[[nm]]], check.keys = TRUE))
    data.frame(construct = nm, n_items = length(constructs[[nm]]), cronbach_alpha = round(a$total$raw_alpha, 3))
}))
cat("\n---- Table S6: Cronbach's alpha ----\n"); print(alpha_tab)

# ---- Tables S5, S7: CFA loadings, AVE, composite reliability ----
cfa_model_si <- '
    PE     =~ PE1 + PE2 + PE3 + PE4
    CP     =~ CP1 + CP2 + CP3 + CP4
    OP     =~ OP1 + OP2 + OP3 + OP4
    BH     =~ BH1 + BH2 + BH3 + BH4
    MT     =~ MT1 + MT2 + MT3 + MT4
    IA_neg =~ IA1 + IA2
    IA_pos =~ IA3 + IA4
    PB     =~ PB1 + PB2 + PB4
    KN     =~ KN1 + KN2 + KN3 + KN4
    TR     =~ TR1 + TR2 + TR3 + TR4
'
fit_cfa_si <- cfa(cfa_model_si, data = data, std.lv = TRUE, estimator = "MLR")
cat("\nCFA (SI) converged:", lavInspect(fit_cfa_si, "converged"), "\n")

kmo_res <- psych::KMO(data[, item_cols])
cat("\n---- KMO measure of sampling adequacy ----\n"); print(round(kmo_res$MSA, 3))
bartlett_res <- psych::cortest.bartlett(cor(data[, item_cols], use = "pairwise.complete.obs"), n = nrow(data))
cat("\n---- Bartlett's test of sphericity ----\n"); print(bartlett_res)

lv_cor <- lavInspect(fit_cfa_si, "cor.lv")
cat("\n---- Latent construct correlation matrix ----\n"); print(round(lv_cor, 3))

load_tab <- standardizedSolution(fit_cfa_si) %>% dplyr::filter(op == "=~") %>%
    dplyr::transmute(construct = lhs, item = rhs, loading = round(est.std,3), pvalue = round(pvalue,4))
cat("\n---- Table S5: standardized CFA loadings ----\n"); print(load_tab)

ave_cr_tab <- load_tab %>% dplyr::group_by(construct) %>%
    dplyr::summarise(n_items = dplyr::n(), AVE = round(mean(loading^2),3),
                     CR = round(sum(loading)^2 / (sum(loading)^2 + sum(1-loading^2)), 3), .groups = "drop")
cat("\n---- Table S6: AVE and composite reliability ----\n"); print(ave_cr_tab)

# ---- Table S8 (measurement-model fit) ----
cfa_fit_indices <- fitMeasures(fit_cfa_si, c("chisq.scaled","df.scaled","pvalue.scaled",
                                             "cfi.robust","rmsea.robust","srmr"))
cat("\n---- Table S8: CFA (measurement model) fit indices ----\n"); print(round(cfa_fit_indices, 4))

# ---- Extended Data 3 (full structural path coefficients, main model) ----
path_tab <- standardizedSolution(sem.fit, type = "std.all", ci = TRUE) %>% dplyr::filter(op == "~") %>%
    dplyr::transmute(lhs, rhs, label, std.est = round(est.std,3), se = round(se,3),
                     ci.lower = round(ci.lower,3), ci.upper = round(ci.upper,3), pvalue = round(pvalue,4))
cat("\n---- Extended Data 3: full structural path coefficients (main model) ----\n"); print(as.data.frame(path_tab))
