# this script is used to compare the results of WMHV models with and without ICV as control variable

library(tidyverse) # version 2.0.0
library(gratia) # version 0.10.0
library(patchwork) # version 1.3.0
library(MetBrewer) # version 0.2.0

path <- "/data/pt_life/ResearchProjects/LLammer/gamms/Results/weighted/"
# load the results of the 1st imputation
load(paste0(path, "imp1/results.RData"))
#define linecol
linecol <- met.brewer(name = "Java", n = 1, type = "discrete")

# get estimates for models 1 and 2 with and without icv as control variable
# and bring them back to the original scale
smicv2 <- smooth_estimates(wmhvsmoothmaxicv) |>
  add_confint(type ="confint")
smicv2 <- transform_fun(smicv2, fun = "exp")
smnoicv2 <- smooth_estimates(wmhvsmoothmax) |>
  add_confint(type ="confint")
smnoicv2 <- transform_fun(smnoicv2, fun = "exp")
ploticv2 <- smicv2 |>
  filter(.smooth == "s(LSNS_SUM)") |> # select the LSNS smooth 
  ggplot() + # initiate plot
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = LSNS_SUM), alpha = 0.1) + # add confint
  geom_line(aes(x = LSNS_SUM, y = .estimate), lwd = 1, color = linecol) + # add smooth line
  labs(y = "WMHV", title = "s(LSNS_SUM) controlled for ICV M2", tag = LETTERS[2]) + theme_bw() +
  scale_x_continuous(limits = c(0, 30), expand = c(0, 0)) 
plotnoicv2 <- smnoicv2 |>
  filter(.smooth == "s(LSNS_SUM)") |> # select the LSNS smooth 
  ggplot() + # initiate plot
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = LSNS_SUM), alpha = 0.1) + # add confint
  geom_line(aes(x = LSNS_SUM, y = .estimate), lwd = 1, color = linecol) + # add smooth line
  labs(y = NULL, title = "s(LSNS_SUM) not controlled for ICV M2", tag = LETTERS[4]) + theme_bw() +
  scale_x_continuous(limits = c(0, 30), expand = c(0, 0)) 
all_plots <- ploticv2 + plotnoicv2
ggsave(paste0(path, "wmhv_with_without_icv.tiff"), all_plots, width = 18, height = 5, units = "in", dpi = 600)

