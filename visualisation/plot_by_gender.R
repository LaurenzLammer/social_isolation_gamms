# this script produces plots depicting the different partial effects of LSNS by gender

library(tidyverse) # version 2.0.0
library(gratia) # version 0.10.0
library(MetBrewer) # version 0.2.0

# define the path
path = "/data/pt_life/ResearchProjects/LLammer/gamms/Results/by_gender/"

#choose palette
palette <- met.brewer(name = "Archambault")

# load the models for each imputation
load(paste0(path, "imp1/models.RData"))

# create one for few and one for many covariates models
matches_min <- list("HCV" = 1, "LESIONS" = 5, "EX. FUNCT." = 9, "MEMORY" = 11, "PROC. SPEED" = 13, 
                    "GAD7_SUM" = 15, "CES_D_SUM" = 17)
matches_max <- list("HCV" = 2, "LESIONS" = 6, "EX. FUNCT." = 10, "MEMORY" = 12, "PROC. SPEED" = 14, 
                    "GAD7_SUM" = 16, "CES_D_SUM" = 18)

composite_plotting <- function(outcome, link = "identity", row = 1, covariates = "min") {
  # this function creates a custom-made gratia plot of the outcome over the different imputations
  # argument link indicates whether results have to be retransformed into the original scale
  # row indicates which row in an assembled line the plot should be in and tag it accordingly
  # argument covariates determines whether models with few or all covariates should be used
  # calculate smooth estimates and confidence intervals for each imputation and bind them into 1 df
  if (covariates == "min") {
    index <- as.numeric(matches_min[outcome])
    abcd <- LETTERS[(row-1)*2 + 1]
    title <- "s(LSNS_SUM) by gender M1"
  } else {
    index <- as.numeric(matches_max[outcome])
    abcd <- LETTERS[(row-1)*2 + 2]
    title <- "s(LSNS_SUM) by gender M2"
  }
  sm <- smooth_estimates(modellist[[index]]) |>
    add_confint(type ="confint")
  if (outcome %in% c("MEMORY", "PROC. SPEED")) {
    sm[,c(".estimate", ".lower_ci", ".upper_ci")] <- sm[,c(".estimate", ".lower_ci", ".upper_ci")] * (-1)
  }
  if (link == "log"){
    sm <- transform_fun(sm, fun = "exp")
  }
  # create the customized plot
  # create the plots comparing partial effects of LSNS by gender
  gender_comparison_plot <- sm |> 
    filter(.smooth %in% c("s(LSNS_SUM):GENDER0", "s(LSNS_SUM):GENDER1")) |> # select the LSNS smooth 
    ggplot() + # initiate plot
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = LSNS_SUM, colour = .smooth), alpha = 0.05) + # add confint
    geom_line(aes(x = LSNS_SUM, y = .estimate, colour = .smooth), lwd = 0.6) + # add smooth line
    labs(y = element_blank(), title = title, tag = abcd)  +
    theme(axis.title.x=element_blank()) + theme_bw() +
    scale_color_manual(values = c(palette[3], palette[1]), labels = c("female", "male"), name = "gender") +
    scale_x_continuous(limits = c(0, 30), expand = c(0, 0)) 
  if(covariates == "min"){
    gender_comparison_plot <- gender_comparison_plot +
      guides(colour="none") +
      labs(y = outcome)
  } 
  return(gender_comparison_plot)
}

# create, arrange and save the plots
HCV_min <- composite_plotting(outcome = "HCV")
LESIONS_min <- composite_plotting(outcome = "LESIONS", link = "log", row = 2)
mri_comp_plots_min <- HCV_min / LESIONS_min
ggsave(paste0(path, "mri_comp_plots_min.tiff"), mri_comp_plots_min, width = 6, height = 8.6, units = "in", dpi = 600)
HCV_max <- composite_plotting(outcome = "HCV", covariates = "max")
LESIONS_max <- composite_plotting(outcome = "LESIONS", link = "log", row = 2, covariates = "max")
mri_comp_plots_max <- HCV_max / LESIONS_max
ggsave(paste0(path, "mri_comp_plots_max.tiff"), mri_comp_plots_max, width = 6, height = 8.6, units = "in", dpi = 600)
mri_comp_plots_min_max <- (HCV_min + HCV_max) / (LESIONS_min + LESIONS_max)
ggsave(paste0(path, "mri_comp_plots_min_max.tiff"), mri_comp_plots_min_max, width = 12, height = 8.6, units = "in", dpi = 600)
EF_min <- composite_plotting(outcome = "EX. FUNCT.")
MEMO_min <- composite_plotting(outcome = "MEMORY", link = "log", row = 2)
PS_min <- composite_plotting(outcome = "PROC. SPEED", row = 3)
cog_comp_plots_min <- EF_min / MEMO_min / PS_min
ggsave(paste0(path, "cog_comp_plots_min.tiff"), cog_comp_plots_min, width = 6, height = 8.6, units = "in", dpi = 600)
EF_max <- composite_plotting(outcome = "EX. FUNCT.", covariates = "max")
MEMO_max <- composite_plotting(outcome = "MEMORY", link = "log", row = 2, covariates = "max")
PS_max <- composite_plotting(outcome = "PROC. SPEED", row = 3, covariates = "max")
cog_comp_plots_max <- EF_max / MEMO_max / PS_max
ggsave(paste0(path, "cog_comp_plots_max.tiff"), cog_comp_plots_max, width = 6, height = 8.6, units = "in", dpi = 600)
cog_comp_plots_min_max <- (EF_min + EF_max) / (MEMO_min + MEMO_max) / (PS_min + PS_max)
ggsave(paste0(path, "cog_comp_plots_min_max.tiff"), cog_comp_plots_min_max, width = 12, height = 8.6, units = "in", dpi = 600)
GAD7_min <- composite_plotting(outcome = "GAD7_SUM", link = "log")
CES_D_min <- composite_plotting(outcome = "CES_D_SUM", link = "log", row = 2)
psy_comp_plots_min <- GAD7_min / CES_D_min 
ggsave(paste0(path, "psy_comp_plots_min.tiff"), psy_comp_plots_min, width = 6, height = 8.6, units = "in", dpi = 600)
GAD7_max <- composite_plotting(outcome = "GAD7_SUM", link = "log", covariates = "max")
CES_D_max <- composite_plotting(outcome = "CES_D_SUM", link = "log", row = 2, covariates = "max")
psy_comp_plots_max <- GAD7_max / CES_D_max
ggsave(paste0(path, "psy_comp_plots_max.tiff"), psy_comp_plots_max, width = 6, height = 8.6, units = "in", dpi = 600)
psy_comp_plots_min_max <- (GAD7_min + GAD7_max) / (CES_D_min + CES_D_max)
ggsave(paste0(path, "psy_comp_plots_min_max.tiff"), psy_comp_plots_min_max, width = 12, height = 8.6, units = "in", dpi = 600)
