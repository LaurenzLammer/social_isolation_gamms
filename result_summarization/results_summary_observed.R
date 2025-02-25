# this script's purpose is to summarise the results of our models run on only observed data
library(tidyverse) # version 2.0.0
library(gratia) # version 0.10.0
library(qvalue) # version 2.38.0

"%not_in%" <- Negate("%in%") # define not in operator

path = "/data/pt_life/ResearchProjects/LLammer/gamms/Results/observed_only/"

# load the workspace to access the overall results table
load(paste0(path, "results.RData"))
# turn the p-values into one-sided p-values
# we must check if the estimate is in the predicted direction (more/less social contact has a beneficial/adverse effect)
# we predicted a positive effect of higher LSNS scores (more social contact) on HCV and EF
# we predicted a negative effect of higher LSNS scores (more social contact) on Lesions, GAD7, CES_D
# the predicted effects on Memory and processing speed are positive in principle
# but as these cognitive functions variables were negated for modelling purposes the predicted effects on them are negative
results <- results[1:36,] # only keep relevant rows
results <- results %>%
  group_by(outcome, covariates) %>%
  mutate(`s_lsns_p-value_onesided` = case_when(
    outcome %in% c("HCV_ADJ", "EF") & lead(LSNS_SUM_Estimate) > 0 ~ 0.5 * `s_lsns_p-value`,
    outcome %in% c("HCV_ADJ", "EF") & lead(LSNS_SUM_Estimate) < 0 ~ 1 - 0.5 * `s_lsns_p-value`,
    outcome %in% c("LESIONS", "LESIONSLOG", "LESIONSLOGICV", "GAD7_SUM", "CES_D_SUM", "MEMO", "PS") & lead(LSNS_SUM_Estimate) < 0 ~ 0.5 * `s_lsns_p-value`,
    outcome %in% c("LESIONS", "LESIONSLOG", "LESIONSLOGICV", "GAD7_SUM", "CES_D_SUM", "MEMO", "PS") & lead(LSNS_SUM_Estimate) > 0 ~ 1 - 0.5 * `s_lsns_p-value`,
  ))
# calculate q-values
results$`s_lsns_q-value` <- NA
results[ results$covariates == 1 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_q-value"] <-
    qvalue(p = as.data.frame(results[results$covariates == 1 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_p-value_onesided"]),
           fdr.level = 0.05, pi0 = 1)$qvalues
results[results$covariates == 2 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_q-value"] <-
    qvalue(p = as.data.frame(results[results$covariates == 2 & results$model == "smooth" & results$outcome %not_in% c("LESIONS", "LESIONSLOG", "CES_D_SUM", "GAD7_SUM"), "s_lsns_p-value_onesided"]),
           fdr.level = 0.05, pi0 = 1)$qvalues
# load the csv files the anova results have been saved to
anovaresults <- read.csv(file = paste0(path, "/anovaresults.csv"))
# calculate the difference in AIC and BIC between smooth and the linear model
anovaresults <- anovaresults %>%
  group_by(outcome, covariates) %>%
  mutate(diff_AIC = lag(AIC) - AIC,
         diff_BIC = lag(BIC) - BIC) 
# create an overview df to collect results
overview <- data.frame(matrix(ncol = 13, nrow = 18))
colnames(overview) <- c("outcome", "covariates", "s_age_edf", "s_age_p-value", "s_lsns_edf", "s_lsns_p-value_onesided", 
                        "s_lsns_q-value", "prop_of_effect", "diff_AIC", "diff_BIC", "Pr..Chisq.", "favours_nonlinear", "nobs")
outcomes <- c("HCV_ADJ", "LESIONS", "LESIONSLOG", "LESIONSLOGICV", "MEMO", "EF", "PS", "GAD7_SUM", "CES_D_SUM")
overview$outcome <- rep(outcomes, each = 2, length.out = 9*2)
overview$covariates <- rep(c("min", "max"), each = 1, length.out = 9*2)

# specify function for data retrieval from results table
get_edf_n_p <- function(outcome, measure){
  # this function returns the value of a measure of the first imputation 
  # along with the minimum and maximum value across all imputations
  overview[overview$outcome == outcome & overview$covariates == "min", measure] <<- 
    results[results$outcome == outcome & results$model == "smooth" & results$covariates == 1, measure]
  overview[overview$outcome == outcome & overview$covariates == "max",measure] <<- 
    results[results$outcome == outcome & results$model == "smooth" & results$covariates == 2, measure]
}
get_anova_measures <- function(outcome, measure){
  overview[overview$outcome == outcome & overview$covariates == "min", measure] <<- 
    anovaresults[anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "min", measure]
  overview[overview$outcome == outcome & overview$covariates == "max", measure] <<- 
    anovaresults[anovaresults$outcome == outcome & anovaresults$model == "smooth" & anovaresults$covariates == "max", measure]
}
return_prop_of_effect <- function(){
  relative_effects <- data.frame(matrix(ncol = 9, nrow = 18))
  colnames(relative_effects) <- c("outcome", "covariates")
  outcomes <- c("HCV_ADJ", "LESIONS", "LESIONSLOG", "LESIONSLOGICV", "MEMO", "EF", "PS", "GAD7_SUM", "CES_D_SUM")
  relative_effects$outcome <- rep(outcomes, each = 2, length.out = 9*2)
  relative_effects$covariates <- rep(c("min", "max"), each = 1, length.out = 9*2)
  load(paste0(path, "/results.RData"))
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
    relative_effects[m, "size"] <- colSums(sm[1:12, "prop_effect"])/colSums(sm[, "prop_effect"])
  }
  return(relative_effects[,"size"])
}
for (outcome in outcomes) {
  for (measure in c("s_age_edf", "s_age_p-value", "s_lsns_edf", "s_lsns_p-value_onesided", "s_lsns_q-value")) {
    get_edf_n_p(outcome = outcome, measure = measure)
  }
  for (measure in c("diff_AIC", "diff_BIC", "Pr..Chisq.")) {
    get_anova_measures(outcome = outcome, measure = measure)
  }
}
overview[,c("prop_of_effect", "min_prop_of_effect", "max_prop_of_effect")] <- return_prop_of_effect()
overview$favours_nonlinear <- ifelse(overview$diff_BIC > 0 & overview$Pr..Chisq. < 0.05, 1, 0)
# retrieve the number of observations for each row in the overview table
for (n in 1:18) {
  m <- (n-1)*2+1
  overview[n,"nobs"] <- modellist[[m]]$df.null
}

# write the overview table to a csv file
write_csv(overview, file = paste0(path, "overview.csv"))
