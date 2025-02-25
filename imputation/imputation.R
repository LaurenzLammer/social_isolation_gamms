#import relevant packages
library(mice) # version 3.17.0
library(miceadds) # version 3.17-44
library(lattice) # version 0.22-6
library(ggmice) # version 0.1.0
library(tidyverse) # version 2.0.0
library(MetBrewer) # version 0.2.0
library(patchwork) # version 1.3.0

# read in the assembled data for the study
data <- read.csv("/data/pt_life/ResearchProjects/LLammer/gamms/Data/assembled_data.csv")

# reduce df to only relevant columns
df <- data %>%
  select("ID", "tp", "AGE", "GENDER", "LSNS_SUM", "CES_D_SUM", "GAD7_SUM", "MEMO", "EF", "PS", "HCV_ADJ", "LESIONS",
         "HYPERTENSION", "DIABETES", "BMI", "EDUCATION", "SES", "MARRIED", "EMPLOYED", "LIVE_ALONE",  matches("^LSNS_[1-6]$"), matches("^CES_D_([1-9]|1[0-9]|20)$"), matches("^GAD7_[1-7]$"), 
         "RADIO_PROBLEMS", "WMHV_PROBLEMS", "COG_COHORT", "MRI_COHORT", "VENTRICLE_EXPANSION", "SES_WEIGHT", "sbTIV", "HEIGHT", "INCOME_SCORE")


# reduce further for plotting of relevant variables
df_plot <- data %>%
  select("LSNS_SUM", "CES_D_SUM", "GAD7_SUM", "MEMO", "EF", "PS", "HCV_ADJ", "LESIONS",
         "HYPERTENSION", "DIABETES", "BMI", "EDUCATION", "SES")
# produce a missingness-pattern plot
pattern_plot <- plot_pattern(df_plot, rotate = T, npat = 20, square = F)
ggsave(plot = pattern_plot, 
       filename = "/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/missingness_patterns.tiff", 
       device = "tiff", dpi = 600, width = 8.2, height = 6.85, bg = "white")

# create categorical LSNS variable based on standard cut-off
df$LSNS_CAT <- ifelse(df$LSNS_SUM < 12, 1, 0)
df$LSNS_CAT <- ifelse(is.na(df$LSNS_SUM), NA, df$LSNS_CAT) # keep NAs
df$LSNS_CAT <- as.factor(df$LSNS_CAT)

# intracranial volume should not differ between observations
# we can thus impute missing observations from one timepoint by observations from the other timepoint if available
df <- df %>%
  group_by(ID) %>%
  arrange(ID, tp == "fu", tp) %>%
  mutate(
    sbTIV = case_when(
      !is.na(sbTIV) ~ sbTIV,
      tp == "bl" & "bl" %in% tp & "fu" %in% tp ~ lead(sbTIV),
      tp == "fu" & "bl" %in% tp & "fu" %in% tp ~ lag(sbTIV),
      TRUE ~ NA_real_  # Otherwise, leave it as NA if the other timepoint is missing
      )) %>%
  ungroup()

predMat <- make.predictorMatrix(df)
# 0 indicates that a variables isn't used for prediction
# 1 indicates that a variable is used as a fixed effect for prediction
# -2 indicates the cluster variable
predMat[c("ID", "tp", "GENDER", "RADIO_PROBLEMS", "WMHV_PROBLEMS", "COG_COHORT", "MRI_COHORT", "VENTRICLE_EXPANSION", "SES_WEIGHT", "HEIGHT"),] <- 0 # define variables that will not be imputed
# let's now prepare the imputation of the questionnaires items
predMat[grepl("^LSNS_[1-6]$", rownames(predMat)) | grepl("^GAD7_[1-7]$", rownames(predMat)) | grepl("^CES_D_\\d+$", rownames(predMat)), ] <- 0

LSNS_cols_to_predict <- c(paste0("LSNS_", c(1:6)))
LSNS_predictors <- c("CES_D_SUM", "AGE", "GENDER", "HYPERTENSION", "BMI", "DIABETES", 
                     "MARRIED", "EMPLOYED", "SES", "LIVE_ALONE")

# Set up the predictor matrix such that LSNS_1 through LSNS_6 are predicted by:
# - Other variables (CESD_SUM, AGE, etc.)
# Set up prediction for LSNS_1 through LSNS_6 to be predicted by the external variables (CESD_SUM, AGE, etc.)
for (col in LSNS_cols_to_predict) {
  predMat[col, LSNS_predictors] <- 1  # These columns will be predicted by the external variables
  predMat[col, "ID"] <- -2 # set ID as grouping variable
}
#repeat for GAD7
GAD7_cols_to_predict <- c(paste0("GAD7_", c(1:7)))
GAD7_predictors <- c("LSNS_SUM", "LSNS_CAT", "AGE", "GENDER", "HYPERTENSION", "BMI", "DIABETES", 
                     "SES", "MARRIED", "EMPLOYED")
for (col in GAD7_cols_to_predict) {
  predMat[col, GAD7_predictors] <- 1  # These columns will be predicted by the external variables
  predMat[col, 1] <- -2 # set ID as grouping variable
}
#repeat for CES_D
CES_D_cols_to_predict <- c(paste0("CES_D_", c(1:20)))
CES_D_predictors <- c("LSNS_SUM", "LSNS_CAT", "AGE", "GENDER", "HYPERTENSION", "BMI", "DIABETES", 
                      "SES", "EMPLOYED")
for (col in CES_D_cols_to_predict) {
  predMat[col, CES_D_predictors] <- 1  # These columns will be predicted by the external variables
  predMat[col, 1] <- -2 # set ID as grouping variable
}
# define predictor matrix
predMat["LSNS_SUM",] <- c(-2, rep(0, 62))
predMat["LSNS_SUM", c(LSNS_predictors, paste0("LSNS_", c(1:6)))] <- 1
predMat["LSNS_CAT",] <- c(-2, rep(0, 62))
predMat["LSNS_CAT", c(LSNS_predictors, paste0("LSNS_", c(1:6)))] <- 1
predMat["CES_D_SUM",] <- c(-2, rep(0, 62))
predMat["CES_D_SUM", c(CES_D_predictors, paste0("CES_D_", c(1:20)))] <- 1
predMat["GAD7_SUM",] <- c(-2, rep(0, 62))
predMat["GAD7_SUM", c(GAD7_predictors, paste0("GAD7_", c(1:7)))] <- 1
predMat["MEMO",] <- c(-2, 0, rep(1, 3), rep(0, 3), rep(1, 8), rep(0, 46), 1)
predMat["EF",] <- c(-2, 0, rep(1, 3), 0, 0, 1, 0, rep(1, 7), rep(0,46), 1)
predMat["PS",] <- c(-2, 0, rep(1, 3), 0, 0, 1, 1, 0, rep(1, 6), rep(0,46), 1)
predMat["HCV_ADJ",] <- c(-2, 0, rep(1, 3), 0, 0, rep(1, 3), 0, rep(1, 5), rep(0,46), 1)
predMat["LESIONS",] <- c(-2, 0, rep(1, 3), 0, 0, rep(1, 4), 0, rep(1, 4), rep(0,46), 1)
predMat["HYPERTENSION",] <- c(0, 0, rep(1, 3), 0, 0, rep(0, 5), rep(0,4), 1, rep(0,46)) # imputed without the grouping variable as BL mostly equals FU, imputed without neurocognitive predictors
predMat["BMI",] <- c(0, 0, rep(1, 2), rep(0, 12), 1, rep(0,46)) # imputed without the grouping variable as BL mostly equals FU, imputed without neurocognitive and psychosocial predictors
predMat["DIABETES",] <- c(0, 0, rep(1, 2), rep(0, 12), 1, rep(0,46)) # imputed without the grouping variable as BL mostly equals FU, imputed without neurocognitive and psychosocial predictors
predMat["EDUCATION",] <- c(0, 0, rep(1, 3), 0, 0, rep(1, 5), 0, 1, 1, 0, 1, 0, 1, rep(0,42), 1, 0) # the last two variables must be imputed without the grouping variable as BL equals FU
predMat["SES",] <- c(0, 0, rep(1, 10), 0, rep(1, 3), 0, 0, 1, rep(0,42), 1, 1)
predMat["AGE",] <- c(-2, 1, 0, rep(0, 9), 1, 1, 1, 0, 0, 0, 1, 1, rep(0, 43))
predMat["MARRIED",] <- c(0, 0, rep(1, 2), rep(0, 12), 1, 0, 0, 1, rep(0, 43))
predMat["EMPLOYED",] <- c(0, rep(1, 3), rep(0, 11), 1, 1, rep(0, 44), 1, 0)
predMat["LIVE_ALONE",] <- c(0, rep(1, 3), rep(0, 12), 1, 1, rep(0, 45))
predMat["INCOME_SCORE",] <- c(0, 0, 1, 1, rep(0, 14), 1, rep(0, 44))
predMat["sbTIV",] <- c(0, 0, 1, rep(0, 57), 1, 0, 0)

# define the imputation method
meth <- c(rep("", 2), "2l.pmm", "", rep("2l.pmm", 8), rep("logreg", 2), rep("pmm", 3), rep("logreg", 3),
          rep("2l.pmm", 33), rep("", 6), "pmm", "", "pmm", "2l.bin")

# save a predictor matrix with a column for the imputation method for visualisation
predmatforviz <- as.data.frame(predMat)
predmatforviz[,"meth"] <- meth
write.csv(predmatforviz, "/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/predmatforviz.csv", row.names = T)

# ensure that categorical variables are recognised as such
df$HYPERTENSION <- as.factor(df$HYPERTENSION)
df$DIABETES <- as.factor(df$DIABETES)
df$MARRIED <- as.factor(df$MARRIED)
df$EMPLOYED <- as.factor(df$EMPLOYED)
df$LIVE_ALONE <- as.factor(df$LIVE_ALONE)

# create a "where"-matrix to avoid imputing observations that have been deemed low quality
where_matrix <- make.where(df)
# do not impute missing values that have been turned into NAs due to problems detected in the QA
where_matrix[df$RADIO_PROBLEMS == 1 | df$WMHV_PROBLEMS == 1, "LESIONS"] <- FALSE
where_matrix[df$RADIO_PROBLEMS == 1, "HCV_ADJ"] <- FALSE
# do not impute MRI values if participants where not part of the MRI cohort
where_matrix[df$MRI_COHORT == 0, c("LESIONS", "HCV_ADJ")] <- FALSE
# do not impute cognitive functions if participants where not part of the cognition cohort
where_matrix[df$COG_COHORT == 0, c("MEMO", "EF", "PS")] <- FALSE

# impute 5 datasets
imp <- mice::mice(data = df, m = 5, seed = 1848, maxit = 25, method = meth, where = where_matrix,
                   predictorMatrix = predMat, donors = 25)

# produce diagnostic plots
# use the colourblind-friendly Hokusai3 palette by MetBrewer
palette <- met.brewer("Hokusai3", n = 5)
# create a trace line plot of means and standard deviations
tiff("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/trace_plots.tiff",
     width = 8.2, height = 6.85, units = "in", res = 600)
plot(imp, y=c("LSNS_SUM","GAD7_SUM", "CES_D_SUM", "MEMO","EF", "PS", "HCV_ADJ","LESIONS"), layout=c(2,8), col = palette)
dev.off()
# create stripplots with boxplots to compare observed and imputed data
# define the variables to depict in the 1st plot
variables1 <- c("LSNS_SUM","GAD7_SUM", "CES_D_SUM", "LESIONS")
plots <- list() # prepare empty variable
# Loop over variables and create individual plots
for (var in variables1) {
  p <- ggmice(imp, aes(x = .imp, y = !!sym(var))) +  
    geom_jitter(height = 0, width = 0.25) +
    geom_boxplot(width = 0.5, size = 1, alpha = 0.75, outlier.shape = NA) +
    labs(x = "Imputation number", y = var)  # Adjust label for each variable
  # Append the plot to the list
  plots[[var]] <- p
}
# Combine all plots into one figure using patchwork
combined_plot <- wrap_plots(plots, ncol = 2)  # Arrange plots in 2 columns 

# save the combined plot
ggsave("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/combined_imputed_variables_plot_quest_and_lesions.tiff", combined_plot, width = 10, height = 12, units = "in", dpi = 600)

# repeat for other variables
variables2 <- c("MEMO","EF", "PS", "HCV_ADJ")
plots <- list() # prepare empty variable
# Loop over variables and create individual plots
for (var in variables2) {
  p <- ggmice(imp, aes(x = .imp, y = !!sym(var))) +  
    geom_jitter(height = 0, width = 0.25) +
    geom_boxplot(width = 0.5, size = 1, alpha = 0.75, outlier.shape = NA) +
    labs(x = "Imputation number", y = var)  # Adjust label for each variable
  # Append the plot to the list
  plots[[var]] <- p
}
# Combine all plots into one figure using patchwork
combined_plot <- wrap_plots(plots, ncol = 2)  # Arrange plots in 2 columns 
# save the combined plot
ggsave("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/combined_imputed_variables_plot_neurocog.tiff", combined_plot, width = 10, height = 12, units = "in", dpi = 600)

# save workspace and mids object for later loading
save.image("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/Workspace.RData")
save(imp, file = "/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/Imputations.RData")
