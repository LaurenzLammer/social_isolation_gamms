# this script's purpose is to summarise the results of our models
library(tidyverse) # version 2.0.0
library(gratia) # version 0.10.0
library(qvalue) # version 2.38.0

"%not_in%" <- Negate("%in%") # define not in operator

# define whether SES-weighted or unweighted results should be used
weighted <- 1 # set to 0 for unweighted models
if (weighted == 1) {
  path = "/data/pt_life/ResearchProjects/LLammer/gamms/Results/weighted/"
} else {
  path = "/data/pt_life/ResearchProjects/LLammer/gamms/Results/unweighted/"
}


# load the workspace to access the overall results table
load(paste0(path, "results.RData"))
# turn the p-values into one-sided p-values
# we must check if the estimate is in the predicted direction (more/less social contact has a beneficial/adverse effect)
# we predicted a positive effect of higher LSNS scores (more social contact) on HCV and EF
# we predicted a negative effect of higher LSNS scores (more social contact) on Lesions, GAD7, CES_D
# the predicted effects on Memory and processing speed are positive in principle
# but as these cognitive functions variables were negated for modelling purposes the predicted effects on them are negative
results <- results %>%
  group_by(outcome, imp, covariates) %>%
  mutate(`s_lsns_p-value_onesided` = case_when(
    outcome %in% c("HCV_ADJ", "EF") & lead(LSNS_SUM_Estimate) > 0 ~ 0.5 * `s_lsns_p-value`,
    outcome %in% c("HCV_ADJ", "EF") & lead(LSNS_SUM_Estimate) < 0 ~ 1 - 0.5 * `s_lsns_p-value`,
    outcome %in% c("LESIONS", "LESIONSLOG", "LESIONSLOGICV", "GAD7_SUM", "CES_D_SUM", "MEMO", "PS") & lead(LSNS_SUM_Estimate) < 0 ~ 0.5 * `s_lsns_p-value`,
    outcome %in% c("LESIONS", "LESIONSLOG", "LESIONSLOGICV", "GAD7_SUM", "CES_D_SUM", "MEMO", "PS") & lead(LSNS_SUM_Estimate) > 0 ~ 1 - 0.5 * `s_lsns_p-value`,
  ))
# calculate q-values
results$`s_lsns_q-value` <- NA
for (i in 1:5) {
  results[results$imp == i & results$covariates == 1 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_q-value"] <-
    qvalue(p = as.data.frame(results[results$imp == i & results$covariates == 1 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_p-value_onesided"]),
           fdr.level = 0.05, pi0 = 1)$qvalues
  results[results$imp == i & results$covariates == 2 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_q-value"] <-
    qvalue(p = as.data.frame(results[results$imp == i & results$covariates == 2 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_p-value_onesided"]),
           fdr.level = 0.05, pi0 = 1)$qvalues
}
# load the csv files the anova results have been saved to
anovaresults <- data.frame()
for (n in 1:5){
  anovares <- read.csv(file = paste0(path, "imp", n, "/anovaresults.csv"))
  anovares$imp <- n
  anovaresults <- rbind(anovaresults, anovares)
}
# calculate the difference in AIC and BIC between smooth and the linear model
anovaresults <- anovaresults %>%
  group_by(outcome, imp, covariates) %>%
  mutate(diff_AIC = lag(AIC) - AIC,
         diff_BIC = lag(BIC) - BIC) 
# create an overview df to collect results
overview <- data.frame(matrix(ncol = 39, nrow = 18))
colnames(overview) <- c("outcome", "covariates", "s_age_edf", "min_s_age_edf", 
                        "max_s_age_edf", "s_age_p-value", "min_s_age_p-value", 
                        "max_s_age_p-value", "s_lsns_edf", "min_s_lsns_edf", 
                        "max_s_lsns_edf", "s_lsns_F", "min_s_lsns_F", 
                        "max_s_lsns_F", "s_lsns_p-value_onesided", "min_s_lsns_p-value_onesided", 
                        "max_s_lsns_p-value_onesided", "prop_lsns_p_significant", "s_lsns_q-value", "min_s_lsns_q-value",  
                        "max_s_lsns_q-value", "prop_lsns_q_significant", "prop_of_effect", "min_prop_of_effect", 
                        "max_prop_of_effect", "diff_AIC", "min_diff_AIC", "max_diff_AIC", 
                        "diff_BIC", "min_diff_BIC", "max_diff_BIC", "Pr..Chisq.", 
                        "min_Pr..Chisq.", "max_Pr..Chisq.", "Chisq", "min_Chisq", "max_Chisq", "favours_nonlinear", "nobs")
outcomes <- c("HCV_ADJ", "LESIONS", "LESIONSLOG", "LESIONSLOGICV", "MEMO", "EF", "PS", "GAD7_SUM", "CES_D_SUM")
overview$outcome <- rep(outcomes, each = 2, length.out = 9*2)
overview$covariates <- rep(c("min", "max"), each = 1, length.out = 9*2)

# specify function for data retrieval from results table
get_imp1_edf_n_p_with_min_and_max <- function(outcome, measure){
  # this function returns the value of a measure of the first imputation 
  # along with the minimum and maximum value across all imputations
  overview[overview$outcome == outcome & overview$covariates == "min", measure] <<- 
    results[results$imp == 1 & results$outcome == outcome & results$model == "smooth" & results$covariates == 1, measure]
  overview[overview$outcome == outcome & overview$covariates == "min", paste0("min_", measure)] <<- 
    min(results[results$outcome == outcome & results$model == "smooth" & results$covariates == 1, measure])
  overview[overview$outcome == outcome & overview$covariates == "min", paste0("max_", measure)] <<- 
    max(results[results$outcome == outcome & results$model == "smooth" & results$covariates == 1, measure])
  overview[overview$outcome == outcome & overview$covariates == "max",measure] <<- 
    results[results$imp == 1 & results$outcome == outcome & results$model == "smooth" & results$covariates == 2, measure]
  overview[overview$outcome == outcome & overview$covariates == "max", paste0("min_", measure)] <<- 
    min(results[results$outcome == outcome & results$model == "smooth" & results$covariates == 2, measure])
  overview[overview$outcome == outcome & overview$covariates == "max", paste0("max_", measure)] <<- 
    max(results[results$outcome == outcome & results$model == "smooth" & results$covariates == 2, measure])
}
get_imp1_anova_measures_with_min_and_max <- function(outcome, measure){
  overview[overview$outcome == outcome & overview$covariates == "min", measure] <<- 
    anovaresults[anovaresults$imp == 1 & anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "min", measure]
  overview[overview$outcome == outcome & overview$covariates == "min", paste0("min_", measure)] <<- 
    min(anovaresults[anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "min", measure])
  overview[overview$outcome == outcome & overview$covariates == "min", paste0("max_", measure)] <<- 
    max(anovaresults[anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "min", measure])
  overview[overview$outcome == outcome & overview$covariates == "max", measure] <<- 
    anovaresults[anovaresults$imp == 1 & anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "max", measure]
  overview[overview$outcome == outcome & overview$covariates == "max", paste0("min_", measure)] <<- 
    min(anovaresults[anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "max", measure])
  overview[overview$outcome == outcome & overview$covariates == "max", paste0("max_", measure)] <<- 
    max(anovaresults[anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "max", measure])
}
return_imp1_prop_of_effect_with_min_max <- function(){
  relative_effects <- data.frame(matrix(ncol = 9, nrow = 18))
  colnames(relative_effects) <- c("outcome", "covariates", paste0("imp", 1:5), "min", "max")
  outcomes <- c("HCV_ADJ", "LESIONS", "LESIONSLOG", "LESIONSLOGICV", "MEMO", "EF", "PS", "GAD7_SUM", "CES_D_SUM")
  relative_effects$outcome <- rep(outcomes, each = 2, length.out = 9*2)
  relative_effects$covariates <- rep(c("min", "max"), each = 1, length.out = 9*2)
  for (n in 1:5) {
    load(paste0(path, "imp", n, "/results.RData"))
    count_table <- table(df$LSNS_SUM)
    proportion_table <- count_table / sum(count_table)
    modellist <- list(hcvsmoothmin$gam, hcvsmoothmax$gam,
                      wmhvsmoothminprereg$gam, wmhvsmoothmaxprereg$gam,
                      wmhvsmoothmin$gam, wmhvsmoothmax$gam,
                      wmhvsmoothminicv$gam, wmhvsmoothmaxicv$gam,
                      EFsmoothmin$gam, EFsmoothmax$gam,
                      MEMOsmoothmin$gam, MEMOsmoothmax$gam,
                      PSsmoothmin$gam, PSsmoothmax$gam,
                      GAD7_SUMsmoothmin$gam, GAD7_SUMsmoothmax$gam,
                      CES_D_SUMsmoothmin$gam, CES_D_SUMsmoothmax$gam)
    for (m in 1:18) {
      sm <- smooth_estimates(modellist[[m]], n = 31)
      sm <- sm[32:62,]
      if (m %in% c(3:8,11,12,15:18)) { # transform estimates if a log transformation has been applied
        sm <- transform_fun(sm, fun = "exp")
      }
      sm$proportion <- proportion_table
      sm <- sm %>%
        mutate(gradient = (abs(.estimate - lag(.estimate, default = first(.estimate))) + abs(.estimate - lead(.estimate, default = last(.estimate)))) / 2)
      sm[c(1,31), "gradient"] <- sm[c(1,31), "gradient"]*2
      sm$prop_effect <- sm$gradient*sm$proportion
      relative_effects[m, paste0("imp", n)] <- colSums(sm[1:12, "prop_effect"])/colSums(sm[, "prop_effect"])
    }
  }
  relative_effects$min <- apply(relative_effects[,3:7], 1, min)
  relative_effects$max <- apply(relative_effects[,3:7], 1, max)
  return(relative_effects[,c("imp1", "min", "max")])
}
for (outcome in outcomes) {
  for (measure in c("s_age_edf", "s_age_p-value", "s_lsns_edf", "s_lsns_p-value_onesided", "s_lsns_q-value",
                    "s_lsns_F")) {
    get_imp1_edf_n_p_with_min_and_max(outcome = outcome, measure = measure)
  }
  for (measure in c("diff_AIC", "diff_BIC", "Pr..Chisq.", "Chisq")) {
    get_imp1_anova_measures_with_min_and_max(outcome = outcome, measure = measure)
  }
  overview[overview$outcome == outcome & overview$covariates == "min", "prop_lsns_p_significant"] <- 
    nrow(results[results$outcome == outcome & results$covariates == 1 & results$model == "smooth" & results$`s_lsns_p-value_onesided` < 0.05,]) /
    nrow(results[results$outcome == outcome & results$covariates == 1 & results$model == "smooth",])
  overview[overview$outcome == outcome & overview$covariates == "max", "prop_lsns_p_significant"] <- 
    nrow(results[results$outcome == outcome & results$covariates == 2 & results$model == "smooth" & results$`s_lsns_p-value_onesided` < 0.05,]) /
    nrow(results[results$outcome == outcome & results$covariates == 2 & results$model == "smooth",])
  overview[overview$outcome == outcome & overview$covariates == "min", "prop_lsns_q_significant"] <- 
    nrow(results[results$outcome == outcome & results$covariates == 1 & results$model == "smooth" & results$`s_lsns_q-value` < 0.05,]) /
    nrow(results[results$outcome == outcome & results$covariates == 1 & results$model == "smooth",])
  overview[overview$outcome == outcome & overview$covariates == "max", "prop_lsns_q_significant"] <- 
    nrow(results[results$outcome == outcome & results$covariates == 2 & results$model == "smooth" & results$`s_lsns_q-value` < 0.05,]) /
    nrow(results[results$outcome == outcome & results$covariates == 2 & results$model == "smooth",])
}
overview[,c("prop_of_effect", "min_prop_of_effect", "max_prop_of_effect")] <- return_imp1_prop_of_effect_with_min_max()
overview$favours_nonlinear <- ifelse(overview$diff_BIC > 0 & overview$Pr..Chisq. < 0.05, 1, 0)
# retrieve the number of observations for each row in the overview table
for (n in 1:18) {
  m <- (n-1)*2+1
  overview[n,"nobs"] <- modellist[[m]]$df.null
}

# write the overview table to a csv file
write_csv(overview, file = paste0(path, "overview.csv"))

# also save the diagnostics table
write_csv(diagnostics, file = paste0(path, "diagnostics.csv"))
