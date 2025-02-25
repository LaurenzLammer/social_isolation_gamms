# this script runs weighted GAMMS on the imputed datasets producing one partial effect of LSNS per gender

# load required packages
library("mgcv") # version 1.9-1
library("mgcv.helper") # version 0.1.9
library("gamm4") # version 0.2-6
library("mice") # version 3.17-44
library("tidyverse") # version 2.0.0

# define whether SES-weighted or unweighted models should be used
path = "/data/pt_life/ResearchProjects/LLammer/gamms/Results/by_gender/"

# prepare a df "diagnostics" to collect all diagnsotic information
diagnostics <- data.frame(matrix(nrow = 5*9*2, ncol = 14)) # nrow = n_imputations * n_outcomes * n_covariate_sets * n_models
colnames(diagnostics) <- c("imp",	"outcome",	"covar",	"k_age",	"edf_age",   "k_index_age",
                           "p_value_age",	"k_lsns", "edf_lsns",   "k_index_lsns",		"p_value_lsns",
                           "concurvity_worst_para",	"concurvity_worst_age",	"concurvity_worst_lsns")
diagnostics$imp <- rep(1:5, each = 36, length.out = 5*9*2)
diagnostics$outcome <- rep(c("HCV_ADJ", "LESIONS", "LESIONSLOG", "LESIONSLOGICV", "MEMO_NEG", "EF", "PS_NEG", "GAD7_SUM", "CES_D_SUM"), each = 4, length.out = 5*9*2)
diagnostics$covar <- rep(c("min", "max"), each = 2, length.out = 5*9*2)

# prepare a df "results" to collect all result information of individual models
results <- data.frame(matrix(nrow = 5*9*2, ncol = 84)) # nrow = n_imputations * n_outcomes * n_covariate_sets * n_models
linear_covariates <- c("(Intercept)", "GENDER", "LSNS_SUM", "AGE", "EDUCATION", "BMI", "DIABETES1", "HYPERTENSION1", "CES_D_SUM_ASINH", "SES", "VENTRICLE_EXPANSION", "sbTIV")
linear_covariate_measures <- c("_Estimate", "_2.50%", "_97.50%", "_Std. Error", "_t_value", "_p")
#colnames(results) <- c("imputation",	"outcome",	"covariates",	"model",	paste0("s_age_", c("edf", "Ref.df", "F", "p-value")),
#                         paste0("s_lsns_", c("edf", "Ref.df", "F", "p-value")), paste0(outer(linear_covariates, linear_covariate_measures, FUN = "paste0")))
colnames(results) <- c("imp",	"outcome",	"covariates",	paste0("s_age_", c("edf", "Ref.df", "F", "p-value")),
                       paste0("s_lsns_", c("edf", "Ref.df", "F", "p-value")), unlist(mapply(function(x, suffix) paste0(x, suffix), rep(linear_covariates, each = length(linear_covariate_measures)), linear_covariate_measures)))
results$imp <- rep(1:5, each = 36, length.out = 5*9*2)
results$outcome <- rep(c("HCV_ADJ", "LESIONS", "LESIONSLOG", "LESIONSLOGICV", "MEMO", "EF", "PS", "GAD7_SUM", "CES_D_SUM"), each = 4, length.out = 5*9*2)
results$covariates <- rep(1:2, each = 2, length.out = 5*9*2)


# specify functions
transform_variables <- function(){
  # a function to transform variables
  # has to be done for every of the imputed datasets
  # turn education into a categorical variable
  df$EDUCATION <- ifelse(df$EDUCATION < 3.6, 1, 0) 
  # LESIONS has to be log-transformed for good model performance
  df$LESIONSLOG <- log(df$LESIONS)
  df$LESIONSLOGICV <- df$LESIONSLOG
  # negate processing speed and memory to turn them from left-skewed to right skewed
  # ensure that there are no values =< 0 so that a gamma distribution can be used
  df$PS_NEG <- df$PS * (-1)
  df$PS_NEG <- df$PS_NEG + 0.00001 - min(df$PS_NEG, na.rm = T)
  df$MEMO_NEG <- df$MEMO * (-1)
  df$MEMO_NEG <- df$MEMO_NEG + 0.00001 - min(df$MEMO_NEG, na.rm = T)
  df$CES_D_SUM_ASINH <- asinh(df$CES_D_SUM)
  df$GENDER <- as.factor(df$GENDER)
  return(df[,c("EDUCATION", "PS_NEG", "MEMO_NEG", "LESIONSLOG", "LESIONSLOGICV", "CES_D_SUM", "GENDER")])
}

formula_creator <- function(outvar, covariates = "min", k_age = 10, k_lsns = 10, icv = 0) {
  # a function to create formulas for GAMMs
  # outvar defines the outcome variable - e.g. HCV_ADJ or EF
  # covariates is either "min" (only age and gender) or "max" (all covariates)
  # k_age determines the k argument for s(AGE)
  # k_lsns determines the k argument for s(LSNS_SUM)
  # icv indicates whether wmhv should be controlled for icv (1 --> control for icv)
  if (k_age == 10 & k_lsns == 10) {
      formula <- paste0(outvar, " ~ s(AGE) + GENDER + s(LSNS_SUM, by = GENDER)")
    } else {
      formula <- paste0(outvar, " ~ s(AGE, k = ", k_age, ") + GENDER + s(LSNS_SUM, k = ", k_lsns, ", by = GENDER)")
    }
  if (covariates == "max") {
    if (outvar %in% c("HCV_ADJ", "MEMO", "EF", "PS")){
      formula <- paste0(formula, " + EDUCATION + BMI + DIABETES + HYPERTENSION + CES_D_SUM_ASINH")
    } else if (outvar %in% c("GAD7_SUM", "CES_D_SUM")){
      formula <- paste0(formula, " + BMI + DIABETES + HYPERTENSION + SES")
    } else { # the case for LESIONS
      formula <- paste0(formula, " + EDUCATION + BMI + DIABETES + HYPERTENSION + CES_D_SUM_ASINH + VENTRICLE_EXPANSION")
    }
    if (icv == 1) {
      formula <- paste0(formula, " + sbTIV")
    }
  }
  return(as.formula(formula))
}

calculate_GAMMs <- function(outvar, covariates = "min", k_age = 10, k_lsns = 10, icv = 0, 
                            cohort = "full", family = "gaussian", link = "identity") {
  # a function to calculate (G)AMMs
  # first six arguments like in formula_creator
  # cohort specifies which subset is to be used - either "full", "MRI" or "COG"
  # family and link define family and link function (surprise)
  form <- formula_creator(outvar = outvar, covariates = covariates, k_age = k_age, k_lsns = k_lsns, icv = icv)
  # define variables dependent on the cohort
  if (cohort == "MRI") {
    dat <<- df[df$MRI_COHORT == 1,]
    Weights <<- df[df$MRI_COHORT == 1,]$SES_WEIGHT
  } else if (cohort == "COG") {
    dat <<- df[df$COG_COHORT == 1,]
    Weights <<- df[df$COG_COHORT == 1,]$SES_WEIGHT
  } else {
    dat <<- df
    Weights <<- df$SES_WEIGHT
  }
  # define variables dependent on the family and link function
  if (family == "gaussian") {
    fam <<- gaussian(link = "identity")
  } else if (family == "poisson") {
    fam <<- poisson(link = "log")
  } else if (family == "gamma" & link == "identity")  {
    fam <<- Gamma(link = "identity") 
  } else if (family == "gamma" & link == "inverse")  {
    fam <<- Gamma(link = "inverse") 
  } else if (family == "gamma" & link == "log")  {
    fam <<- Gamma(link = "log") 
  } 
  # find the row-number of the diagnostics and results df
  cond1 <- diagnostics$imp == n
  cond2 <- diagnostics$outcome == outvar
  cond3 <- diagnostics$covar == covariates
  rown <- which(cond1 & cond2 & cond3)
  # calculate the (G)AMMs
  res <- gamm4(data = dat, formula = form, random =~(1|ID), REML = F, weights = Weights, family = fam)
  # get k-check data and save it in diagnostics df
  kcheck <- k.check(res$gam)
  diagnostics[rown, 4:7] <<- kcheck["s(AGE)",]
  try(diagnostics[rown, 8:11] <<- kcheck["s(LSNS_SUM)",], silent = T) # use try as not every model has smooth LSNS term
  # save an image of the diagnostic plots
  tiff(paste0(path, "imp", n, "/diagnostic_plots/",
              outvar, covariates, ".tiff"),
       width = 8.2, height = 8.2, units = "in", res = 600)
  par(mfrow = c(2, 2))
  gam.check(res$gam)
  dev.off()
  # get concurvity data and save it in diagnostics df
  concurv <- as.data.frame(concurvity(res$gam))
  diagnostics[rown, 13:14] <<- concurv[1, 1:2]
  try(diagnostics[rown, 15] <<- concurv[1, 3], silent = T)
  # get a summary of the result of the (G)AMM
  sres <- summary(res$gam)
  # save results for the smooth terms in the df
  results[rown, 4:7] <<- sres$s.table[1, ]
  try(results[rown, 8:11] <<- sres$s.table[2, ], silent = T)
  # return the (G)AMM
  return(res)
}

# load the imputed datasets
load("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/Imputations.RData")

# iterate over the 5 imputed datasets
for(n in 1:5){
  df <- complete(imp, action = n) # extract the nth dataset
  df[,c("EDUCATION", "PS_NEG", "MEMO_NEG", "LESIONSLOG", "LESIONSLOGICV", "CES_D_SUM_ASINH", "GENDER")] <- transform_variables()
  # run all models with their specific arguments
  hcvsmoothmin <- calculate_GAMMs(outvar = "HCV_ADJ", cohort = "MRI")
  hcvsmoothmax <- calculate_GAMMs(outvar = "HCV_ADJ", covariates = "max", cohort = "MRI")
  print(paste0("HCV", n, " done"))
  print(Sys.time())
  wmhvsmoothmin <- calculate_GAMMs(outvar = "LESIONSLOG", cohort = "MRI")
  wmhvsmoothmax <- calculate_GAMMs(outvar = "LESIONSLOG", covariates = "max", cohort = "MRI")
  print(paste0("wmhv", n, " done"))
  print(Sys.time())
  wmhvsmoothminprereg <- calculate_GAMMs(outvar = "LESIONS", cohort = "MRI", family = "gamma", link = "log", k_age = 20)
  wmhvsmoothmaxprereg <- calculate_GAMMs(outvar = "LESIONS", covariates = "max", cohort = "MRI", family = "gamma", link = "log", k_age = 20)
  print(paste0("wmhv as prereg ", n, " done"))
  print(Sys.time())
  wmhvsmoothminicv <- calculate_GAMMs(outvar = "LESIONSLOGICV", cohort = "MRI", icv = 1)
  wmhvsmoothmaxicv <- calculate_GAMMs(outvar = "LESIONSLOGICV", covariates = "max", cohort = "MRI", icv = 1)
  print(paste0("wmhv controlled for ICV ", n, " done"))
  print(Sys.time())
  MEMOsmoothmin <- calculate_GAMMs(outvar = "MEMO_NEG", cohort = "COG", family = "gamma", link = "log")
  MEMOsmoothmax <- calculate_GAMMs(outvar = "MEMO_NEG", covariates = "max", cohort = "COG", family = "gamma", link = "log")
  print(paste0("memo", n, " done"))
  print(Sys.time())
  EFsmoothmin <- calculate_GAMMs(outvar = "EF", cohort = "COG")
  EFsmoothmax <- calculate_GAMMs(outvar = "EF", covariates = "max", cohort = "COG")
  print(paste0("EF", n, " done"))
  print(Sys.time())
  PSsmoothmin <- calculate_GAMMs(outvar = "PS_NEG", k_age = 20, k_lsns = 20, family = "gamma", link = "identity")
  PSsmoothmax <- calculate_GAMMs(outvar = "PS_NEG", covariates = "max", k_age = 20, k_lsns = 20, family = "gamma", link = "identity")
  print(paste0("PS", n, " done"))
  print(Sys.time())
  GAD7_SUMsmoothmin <- calculate_GAMMs(outvar = "GAD7_SUM", family = "poisson", link = "log")
  GAD7_SUMsmoothmax <- calculate_GAMMs(outvar = "GAD7_SUM", covariates = "max", family = "poisson", link = "log")
  print(paste0("GAD7_SUM", n, " done"))
  print(Sys.time())
  CES_D_SUMsmoothmin <- calculate_GAMMs(outvar = "CES_D_SUM", family = "poisson", link = "log")
  CES_D_SUMsmoothmax <- calculate_GAMMs(outvar = "CES_D_SUM", covariates = "max", family = "poisson", link = "log")
  print(paste0("CES_D_SUM", n, " done"))
  print(Sys.time())
  # save the models for this imputation
  modellist <- list(hcvsmoothmin$gam, hcvsmoothmax$gam,
                    wmhvsmoothmin$gam, wmhvsmoothmax$gam,
                    wmhvsmoothminicv$gam, wmhvsmoothmaxicv$gam,
                    wmhvsmoothminprereg$gam, wmhvsmoothmaxprereg$gam,
                    EFsmoothmin$gam, EFsmoothmax$gam,
                    MEMOsmoothmin$gam, MEMOsmoothmax$gam,
                    PSsmoothmin$gam, PSsmoothmax$gam,
                    GAD7_SUMsmoothmin$gam, GAD7_SUMsmoothmax$gam,
                    CES_D_SUMsmoothmin$gam, CES_D_SUMsmoothmax$gam) 
  save(modellist, file = paste0(path, "imp", n, "/models.RData"))
  # save the workspace for this imputation to an RData file
  save.image(paste0(path, "imp", n, "/results.RData"))
}
write_csv(results, paste0(path, "results.csv"))
save.image(paste0(path, "results.RData"))
