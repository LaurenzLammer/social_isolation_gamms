# this script checks VIFs for the linear terms in the models 
library(car) # version 3.1-3

# prepare df to store VIF results
vifs <- data.frame(matrix(nrow = 5*7, ncol = 3)) # nrow = n_imputations * n_outcomes
colnames(vifs) <- c("imp",	"outcome",	"max_VIF")
vifs$imp <- rep(1:5, each = 7, length.out = 5*7)
vifs$outcome <- rep(c("HCV_ADJ", "LESIONS", "MEMO", "EF", "PS", "GAD7_SUM", 
                      "CES_D_SUM"), each = 1, length.out = 5*7)



#load imputed datasets
load("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Imputation/Imputations.RData")
# iterate over the 5 imputed datasets
for(n in 1:5){
  df <- complete(imp, action = n)
  df$EDUCATION <- ifelse(df$EDUCATION < 3.6, 1, 0) 
  df$CES_D_SUM_ASINH <- asinh(df$CES_D_SUM)
  hcv <- lm(data = df, HCV_ADJ ~ GENDER + BMI + EDUCATION + HYPERTENSION + DIABETES)
  lesions <- lm(data = df, LESIONS ~ GENDER + BMI + EDUCATION + HYPERTENSION + DIABETES + sbTIV + VENTRICLE_EXPANSION)
  memo <- lm(data = df, MEMO ~ GENDER + BMI + EDUCATION + HYPERTENSION + DIABETES)
  ef <- lm(data = df, EF ~ GENDER + BMI + EDUCATION + HYPERTENSION + DIABETES)
  ps <- lm(data = df, PS ~ GENDER + BMI + EDUCATION + HYPERTENSION + DIABETES)
  gad7 <- lm(data = df, GAD7_SUM ~ GENDER + SES + BMI + HYPERTENSION + DIABETES)
  cesd <- lm(data = df, CES_D_SUM ~ GENDER + SES + BMI + HYPERTENSION + DIABETES)
  vifs[(n-1)*7+1,"max_VIF"] <- max(vif(hcv))
  vifs[(n-1)*7+2,"max_VIF"] <- max(vif(lesions))
  vifs[(n-1)*7+3,"max_VIF"] <- max(vif(memo))
  vifs[(n-1)*7+4,"max_VIF"] <- max(vif(ef))
  vifs[(n-1)*7+5,"max_VIF"] <- max(vif(ps))
  vifs[(n-1)*7+6,"max_VIF"] <- max(vif(gad7))
  vifs[(n-1)*7+7,"max_VIF"] <- max(vif(cesd))
}
# save the maximum VIFs
write.csv(vifs, "/data/pt_life/ResearchProjects/LLammer/gamms/Results/vifs.csv")
