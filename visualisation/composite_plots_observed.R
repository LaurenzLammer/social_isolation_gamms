# this script produces plots for analyses that relied on observed data only

library(tidyverse) # version 2.0.0
library(gratia) # version 0.10.0
library(patchwork) # version 1.3.0
library(MetBrewer) # version 0.2.0

# this script produces composite plots of the partial effects of LSNS on our outcomes comparing them to the effect of age

path = "/data/pt_life/ResearchProjects/LLammer/gamms/Results/observed_only/"
# load overview or results df
overview <- read.csv(paste0(path, "overview.csv"))
# rename outcome values in overview to match script
overview[c(1:4,7:18), "outcome"] <- rep(c("HCV", "LESIONSNOLOG", "LESIONS", "EX. FUNCT.", "MEMORY", "PROC. SPEED", 
                                        "GAD7_SUM", "CES_D_SUM"), each = 2)

# choose palette & linecolour
linecol <- met.brewer(name = "Java", n = 1, type = "discrete")
palette_density <- met.brewer(name="OKeeffe1", n=2, type="discrete")

# load uncentered data to plot age in years since birth rather than as difference to mean age in years
uncentered <- read.csv("/data/pt_life/ResearchProjects/LLammer/gamms/Data/assembled_unscaled_data.csv")
mean_age <- mean(uncentered$AGE, na.rm = T)

# load the results
load(paste0(path, "/results.RData"))
# to supplement the partial effect of LSNS we want a two-colored density plot 
# illustrating the proportion of persons above and below the standard cut-off of 12
# create a standard density plot as the basis
temporary_plot <- df %>%  
  ggplot() + 
  geom_density(aes(x = LSNS_SUM))
# get the density values out of the plot
build <- ggplot2::ggplot_build(temporary_plot)
# define a factor by the common LSNS cut-off
df_breaks <- build$data[[1]] %>% 
  mutate(status = case_when(x < 12 ~ 'class1', # standard LSNS-6 cut-off
                            TRUE ~ 'class2'))
# create a two-colored density plot with the color being defined by which side of the cut-off the value is on
density_lsns <- df_breaks %>% 
  ggplot() + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
                   legend.position = "none",     panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                   panel.grid.major = element_blank(), #remove major gridlines
                   panel.grid.minor = element_blank()) +
  geom_area(aes(x = x, y = y, fill = status), alpha=0.4) + scale_fill_manual(values=palette_density) +
  scale_x_continuous(limits = c(0, 30), expand = c(0, 0)) 

# load the models 
load(paste0(path, "/models.RData"))

# create a list matching outcomes to the indices of the models saved in imp1, imp2, ...
# create one for few and one for many covariates models
matches_min <- list("HCV" = 1, "LESIONS" = 9, "EX. FUNCT." = 17, "MEMORY" = 21, "PROC. SPEED" = 25, 
                    "GAD7_SUM" = 29, "CES_D_SUM" = 33)
matches_max <- list("HCV" = 3, "LESIONS" = 11, "EX. FUNCT." = 19, "MEMORY" = 23, "PROC. SPEED" = 27, 
                    "GAD7_SUM" = 31, "CES_D_SUM" = 35)

composite_plotting <- function(outcome, link = "identity", row = 1, covariates = "min") {
  # this function creates a custom-made gratia plot of the outcome over the different imputations
  # argument link indicates whether results have to be retransformed into the original scale
  # row indicates which row in an assembled line the plot should be in and tag it accordingly
  # argument covariates determines whether models with few or all covariates should be used
  # extract significance of LSNS smooth for the 1st imputation
  significance <- data.frame(significance = "")
  significance <- significance %>% mutate(significance = case_when(
    overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_q.value"] < 0.001 ~ "***",
    overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_q.value"] < 0.01 ~ "**",
    overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_q.value"] < 0.05 ~ "*", 
    is.na(overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_q.value"]) & overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_p.value_onesided"]  < 0.001 ~ "***",
    is.na(overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_q.value"]) & overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_p.value_onesided"]  < 0.01 ~ "**",
    is.na(overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_q.value"]) & overview[overview$outcome == outcome & overview$covariates == covariates, "s_lsns_p.value_onesided"]  < 0.05 ~ "*",
    TRUE ~ ''))
  # calculate smooth estimates and confidence intervals for each imputation and bind them into 1 df
  if (covariates == "min") {
    index <- as.numeric(matches_min[outcome])
  } else {
    index <- as.numeric(matches_max[outcome])
  }
  sm <- smooth_estimates(modellist[[index]]) |>
    add_confint(type ="confint")
  sm$AGE <- sm$AGE + mean_age
  if (outcome %in% c("MEMORY", "PROC. SPEED")) {
    sm[,c(".estimate", ".lower_ci", ".upper_ci")] <- sm[,c(".estimate", ".lower_ci", ".upper_ci")] * (-1)
  }
  if (link == "log"){
    sm <- transform_fun(sm, fun = "exp")
  }
  # create the customized plot
    # create plot of partial effect of LSNS 
  plotlsns <- sm |>
    filter(.smooth == "s(LSNS_SUM)") |> # select the LSNS smooth and the 1st imputation
    ggplot() + # initiate plot
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = LSNS_SUM), alpha = 0.1) + # add confint
    geom_line(aes(x = LSNS_SUM, y = .estimate), lwd = 1, color = linecol) + # add smooth line
    labs(y = NULL, title = "s(LSNS_SUM)", tag = LETTERS[(row - 1)*2+2]) + theme_bw() +
    geom_vline(xintercept = 12, color="firebrick1", linetype = "dashed") + # add vertical line at LSNS cut-off
    annotate("text", x = 28, y = Inf, label = significance$significance, hjust = "right", vjust = 1.2, size = 8, color = "black") +
    scale_x_continuous(limits = c(0, 30), expand = c(0, 0)) 
  # create plot of partial effect of age 
  plotage <- sm |> # create a plot for the age smooth
    filter(.smooth == "s(AGE)") |> # select the LSNS smooth and the 1st imputation
    ggplot() + # initiate plot
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = AGE), alpha = 0.1) + # add confint
    geom_line(aes(x = AGE, y = .estimate), lwd = 1, color = linecol) + # add smooth line
    labs(y = paste0("Partial effect on ", outcome), title = "s(AGE)", tag = LETTERS[(row - 1)*2+1]) + theme_bw() +
    scale_x_continuous(limits = c(min(sm$AGE, na.rm = T), max(sm$AGE, na.rm = T)), expand = c(0, 0)) 
  # prepare a patchwork layout to overlay the density plot on the lsns plot
  layout_density <- c(
    area(t = 1, l = 1, b = 5, r = 5),
    area(t = 5, l = 1, b = 5, r = 5)
  )
  # combine the two lsns plots
  combined_plot_lsns <- plotlsns + density_lsns + plot_layout(design = layout_density)
  
  # arrange age and lsns plot next to each other
  combined_plots <- plotage + combined_plot_lsns 
  plots <- list(combined_plots, combined_plot_lsns, plotage, plotlsns)
  return(plots)
}

# create, arrange and save the plots
HCV_min <- composite_plotting(outcome = "HCV")
LESIONS_min <- composite_plotting(outcome = "LESIONS", link = "log", row = 2)
mri_comp_plots_min <- HCV_min[[1]] / LESIONS_min[[1]] 
ggsave(paste0(path, "mri_comp_plots_min.tiff"), mri_comp_plots_min, width = 18, height = 8.6, units = "in", dpi = 600)
HCV_max <- composite_plotting(outcome = "HCV", covariates = "max")
LESIONS_max <- composite_plotting(outcome = "LESIONS", link = "log", row = 2, covariates = "max")
mri_comp_plots_max <- HCV_max[[1]] / LESIONS_max[[1]] 
ggsave(paste0(path, "mri_comp_plots_max.tiff"), mri_comp_plots_max, width = 18, height = 8.6, units = "in", dpi = 600)
EF_min <- composite_plotting(outcome = "EX. FUNCT.")
MEMO_min <- composite_plotting(outcome = "MEMORY", link = "log", row = 2)
PS_min <- composite_plotting(outcome = "PROC. SPEED", row = 3)
cog_comp_plots_min <- EF_min[[1]] / MEMO_min[[1]] / PS_min[[1]]
ggsave(paste0(path, "cog_comp_plots_min.tiff"), cog_comp_plots_min, width = 18, height = 8.6, units = "in", dpi = 600)
EF_max <- composite_plotting(outcome = "EX. FUNCT.", covariates = "max")
MEMO_max <- composite_plotting(outcome = "MEMORY", link = "log", row = 2, covariates = "max")
PS_max <- composite_plotting(outcome = "PROC. SPEED", row = 3, covariates = "max")
cog_comp_plots_max <- EF_max[[1]] / MEMO_max[[1]] / PS_max[[1]]
ggsave(paste0(path, "cog_comp_plots_max.tiff"), cog_comp_plots_max, width = 18, height = 8.6, units = "in", dpi = 600)
GAD7_min <- composite_plotting(outcome = "GAD7_SUM", link = "log")
CES_D_min <- composite_plotting(outcome = "CES_D_SUM", link = "log", row = 2)
psy_comp_plots_min <- GAD7_min[[1]] / CES_D_min[[1]] 
ggsave(paste0(path, "psy_comp_plots_min.tiff"), psy_comp_plots_min, width = 18, height = 8.6, units = "in", dpi = 600)
GAD7_max <- composite_plotting(outcome = "GAD7_SUM", link = "log", covariates = "max")
CES_D_max <- composite_plotting(outcome = "CES_D_SUM", link = "log", row = 2, covariates = "max")
psy_comp_plots_max <- GAD7_max[[1]] / CES_D_max[[1]]
ggsave(paste0(path, "psy_comp_plots_max.tiff"), psy_comp_plots_max, width = 18, height = 8.6, units = "in", dpi = 600)
