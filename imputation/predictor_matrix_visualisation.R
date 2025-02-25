# this script serves to create a visualisation of the imputation predictor matrix and methods
library(tidyverse) # version 2.0.0
library(MetBrewer) # version 0.2.0

# define colourblind-frinedly palette
palet <- met.brewer("Hiroshige")
# load predictormatrix with method column as df
predmat <- read.csv("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/predmatforviz.csv", row.names = 1)
# only keep one row/col per questionnaire for the item imputation and remove irrelevant rows/cols
predmat <- predmat[!rownames(predmat) %in% c(paste0("LSNS_", 2:6), paste0("CES_D_", 2:20), paste0("GAD7_", 2:7), "RADIO_PROBLEMS", "WMHV_PROBLEMS", "COG_COHORT", "MRI_COHORT", "VENTRICLE_EXPANSION", "SES_WEIGHT", "HEIGHT"),
                                   !colnames(predmat) %in% c(paste0("LSNS_", 2:6), paste0("CES_D_", 2:20), paste0("GAD7_", 2:7), "RADIO_PROBLEMS", "WMHV_PROBLEMS", "COG_COHORT", "MRI_COHORT", "VENTRICLE_EXPANSION", "SES_WEIGHT", "HEIGHT")]
# rename rows and cols
rownames(predmat)[rownames(predmat) == "LSNS_1"] <- "LSNS items"
rownames(predmat)[rownames(predmat) == "CES_D_1"] <- "CES_D items"
rownames(predmat)[rownames(predmat) == "GAD7_1"] <- "GAD7 items"
colnames(predmat)[colnames(predmat) == "LSNS_1"] <- "LSNS items"
colnames(predmat)[colnames(predmat) == "CES_D_1"] <- "CES_D items"
colnames(predmat)[colnames(predmat) == "GAD7_1"] <- "GAD7 items"

# get predmat without method col
predmat_no_meth <- predmat[,1:26]

# turn into long format
predmat_melted <- melt(as.matrix(predmat_no_meth))
# name values
predmat_melted$value <- factor(predmat_melted$value, levels = c(-2, 0, 1), labels = c("Grouping" ,"Not Used", "Used"))
# reintegrate imputation method
predmat_melted$Method <- rep(predmat[,27], each = nrow(predmat_no_meth))
# only keep information on method once to avoid cluttering the plot
predmat_melted$Method <- ifelse(predmat_melted$Var1 == "ID", predmat_melted$Method, "")

# visualise matrix with imputation method and save it
viz <- ggplot(predmat_melted, aes(Var2, Var1)) +
  # Predictor matrix as heatmap
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_manual(values = c(palet[1], palet[5], palet[10])) +  
  # Add imputation method as text label or additional layer
  geom_text(aes(label = Method), color = "black", size = 2.5) +
  theme_minimal() +
  labs(
    title = "Predictor Matrix with Imputation Method",
    x = "Variables",
    y = "Variables"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# save the figure
ggsave(plot = viz, filename = "/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/predmatviz.tiff", width = 17, height = 8, units = "in", dpi = 600, device = "tiff")
