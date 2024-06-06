# set working directory
setwd("/data/pt_life/ResearchProjects/LLammer/gamms/Analysis/Simulations")
# load required package
library(gamm4)
# define functions for LME modelling
simulate_amm <- function(dv){ # dv is respective dependent variable
  fullmod <- gamm4(data = df[df$neurocog == 1,], formula = as.formula(paste0(dv, " ~ s(age) + gender + s(LSNS)")), random =~(1|ID), REML = F)
  basemod <- gamm4(data = df[df$neurocog == 1,], formula = as.formula(paste0(dv, " ~ s(age) + gender + LSNS")), random =~(1|ID), REML = F)
  # save plot of s(LSNS) for qualitative inspection
  # Open device
  tiff(paste0("plots/", dv,n), res = 600, width = 8, height = 8, units = "cm")
  # Make a plot
  plot(fullmod$gam, select = 2)
  # Close device
  dev.off()
  # compare full and base model
  anovares <- anova(basemod$mer, fullmod$mer)
  # create a summary of the model to return the edf of the LSNS-smoot
  modelsum <- summary(fullmod$gam)
  # save results into external dataframe
  return(c(anovares[1,c("AIC", "BIC", "Pr(>Chisq)")], anovares[2,c("AIC", "BIC", "Pr(>Chisq)")], modelsum$s.table[2,1]))
}
simulate_gamm <- function(dv){
  fullmod <- gamm4(data = df, formula = as.formula(paste0(dv, " ~ s(age) + gender + s(LSNS)")), random =~(1|ID), family = poisson(link = "log"), REML = F)
  basemod <- gamm4(data = df, formula = as.formula(paste0(dv, " ~ s(age) + gender + LSNS")), random =~(1|ID), family = poisson(link = "log"), REML = F)
  # use predict to plot s(LSNS) on scale of response variable and save plot for qualitative inspection
  # Open device
  tiff(paste0("plots/", dv,n), res = 600, width = 8, height = 8, units = "cm")
  # Make a plot
  # first, create new dataframe for predict function with age and gender set to 0
  newdata = data.frame(LSNS=0:30, age= 0, gender = 0)
  # plot on response variable scale
  plot(predict(fullmod$gam, newdata, type="response"))
  # Close device
  dev.off()
  # compare full and base model
  anovares <- anova(basemod$mer, fullmod$mer)
  # create a summary of the model to return the edf of the LSNS-smoot
  modelsum <- summary(fullmod$gam)
  # save results into external dataframe
  return(c(anovares[1,c("AIC", "BIC", "Pr(>Chisq)")], anovares[2,c("AIC", "BIC", "Pr(>Chisq)")], modelsum$s.table[2,1]))
}

# prepare data frame for simulation results
simres <- data.frame(matrix(nrow = 100*5*7, ncol = 9))
colnames(simres) <- c("simulation", "outcome", "AICbase", "BICbase", "pbase", "AICfull", "BICfull", "pfull", "edf")
simres$simulation <- unlist(lapply(1:100, rep, times = 5*7))
simres$outcome <- as.vector(outer(c("HCV", "exfunct", "memo", "procspeed", "CESD_out", "WMHV", "GAD"),c("", "exp", "lin", "threshexp", "threshlin"), paste, sep = "" ))

# prepare vectors of outcome variables for later use
outvars_gauss <- as.vector(outer(c("HCV", "exfunct", "memo", "procspeed"),c("", "exp", "lin", "threshexp", "threshlin"), paste, sep = "" ))
outvars_poisson <- as.vector(outer(c("CESD_out", "WMHV", "GAD"),c("", "exp", "lin", "threshexp", "threshlin"), paste, sep = "" ))

# prepare dataframe for simulated variables
df <- data.frame(matrix(ncol = 22, nrow = 10000+5500))
colnames(df) <- c("ID", "fu", "neurocog", "LSNS", "CESD", "BMI", "hypertension", "diabetes", "education", "age", "gender", 
                  "SES", "PSQI", "IPAQ", "weight", "HCV_ri", "WMHV_ri", "memo_ri", "procspeed_ri", "exfunct_ri", 
                  "CESD_ri", "GAD_ri")
df$ID[1:10000] <- as.character(1:10000) # 10k participants at baseline
df$ID[10001:15500] <- as.character(1:5500) # 5k returned for follow-up
df$fu[1:10000] <- 0
df$fu[10001:15500] <- 1
df$neurocog[c(1:2500, 10001:11800)] <- 1 # 2500 of baseline and 1800 of follow-up participants were part of the neurocognitive subcohort
df$neurocog[-c(1:2500, 10001:11800)] <- 0

# loop through 100 simulations
for (n in 1:100) {
  print(n)
  set.seed(n) # set seed for reproducibility
  # simulate predictor variables based on prevalence / mean and SD from Lammer et al.
  df$LSNS <- round(rnorm(15500, mean = 14.1, sd = 5.2), digits = 0)
  # define a variable that is 1 for all above threshold observations and 0 for all below threshold
  df$thresh <- ifelse(df$LSNS > 18, 1, 0)
  df$CESD <- rpois(15500, 10) # just as a predictor not as an outcome
  # CESD as predictor is log-transformed
  df$CESD <- log(df$CESD + 0.00001)
  df$age[1:10000] <- rnorm(10000, mean = 0, sd = 7) # CAVE: mean = 0 cause its centered --- the actual sample will likely be somewhat younger but shouldn't matter
  df$age[10001:15500] <- df$age[1:5500] + rnorm(5500, 5.89, 1.97)
  df$gender[1:10000] <- rbinom(n = 10000, prob = 0.54, size = 1)
  df$gender[10001:15500] <- df$gender[1:5500] # we will assume gender to be constant between BL and FU
  df$PSQI <- rnorm(15500, mean = 6, sd = 3)
  df$IPAQ <- rnorm(15500, mean = 5910, sd = 5829)
  # IPAQ and PSQI are z-scored
  df$PSQI <- as.vector(scale(df$PSQI))
  df$IPAQ <- as.vector(scale(df$IPAQ))
  df$hypertension[1:10000] <- rbinom(n = 10000, prob = 0.61, size = 1)
  df$hypertension[10001:15500] <- df$hypertension[1:5500] # hypertensive at BL will be so at FU
  df$hypertension <- ifelse(df$fu == 1 & df$hypertension == 0 & as.numeric(df$ID) %% 10 == 0, 1, 
                            df$hypertension) # an additional 10% of nonhypertensives will turn hypertensive
  df$education[1:10000] <- rbinom(n = 10000, prob = 0.13, size = 1)
  df$education[10001:15500] <- df$education[1:5500] # we will assume education to be constant between BL and FU
  df$BMI <- rnorm(n = 15500, mean = 27.9, sd = 4.2)
  df$diabetes[1:10000] <- rbinom(n = 10000, prob = 0.18, size = 1)
  df$diabetes[10001:15500] <- df$diabetes[1:5500] # diabetic at BL will be so at FU
  df$diabetes <- ifelse(df$fu == 1 & df$diabetes == 0 & as.numeric(df$ID) %% 20 == 0, 1, 
                        df$hypertension) # an additional 5% of nondiabetics will turn diabetic
  df$SES[1:10000] <- sample((1:5), size = 10000, replace = T, prob = c(6.75, 15.36, 22.01, 23.82, 32.07))
  df$SES[10001:15500] <- df$SES[1:5500] # we will assume SES to be constant between BL and FU
  SESweights <- c(0.2 / (length(which(df$SES == 1))/nrow(df)), 0.2 / (length(which(df$SES == 2))/nrow(df)), 
                  0.2 / (length(which(df$SES == 3))/nrow(df)), 0.2 / (length(which(df$SES == 4))/nrow(df)), 
                  0.2 / (length(which(df$SES == 5))/nrow(df)))
  df$weight <- SESweights[df$SES]
  
  
  # simulate outcome variables under the assumption of a linear relationship 
  # use adjusted regression coefficients of models without interaction terms from Lammer et al. (App, Tab 4) and App. Tabs 13 & 14 for the coefficients for IPAQ and PSQI 
  # intercept and random intercept SDs are based on psqi models with max covariates for HCV and cog functions
  # residual SDs are extracted from the same models
  # HCV: intercept = 3828.409, random intercept sd = 309.07, residual sd = 47.4
  # ef: intercept = 0.643381, random intercept sd= 0.6415, residual sd = 0.466906
  # memo: intercept = 0.479806, random intercept sd= 0.5878, residual sd = 0.37214
  # procspeed: intercept = 0.462252, random intercept sd= 0.4993, residual sd = 0.3964443
  
  # create random intercepts for each subject and each outcome
  # the SDs of the random intercepts are drawn from the data of Lammer et al. 
  df$HCV_ri[1:10000] <- rnorm(10000, 0, 309.07)
  df$HCV_ri[10001:15500] <- df$HCV_ri[1:5500]
  df$exfunct_ri[1:10000] <- rnorm(10000, 0, 0.6415)
  df$exfunct_ri[10001:15500] <- df$exfunct_ri[1:5500]
  df$memo_ri[1:10000] <- rnorm(10000, 0, 0.5878)
  df$memo_ri[10001:15500] <- df$memo_ri[1:5500]
  df$procspeed_ri[1:10000] <- rnorm(10000, 0, 0.4993)
  df$procspeed_ri[10001:15500] <- df$procspeed_ri[1:5500]
  # CESD, WMHV and GAD7 were no dependent variables in Lammer et al.
  # to approximate the SDs of the random intercepts for these variables, we will use the average SD of the 4 vars above relative to the variables' general SD
  # the SD is 1 for the cog functions and ca. 420 (411 at BL and 430 at FU) for HCV
  # SD_ri = (0.6415 + 0.5878 + 0.4993 + 309.07 / 420 ) / 4 = 0.6161202
  
  # the SD for CESD is 6 as reported in Lammer et al.
  # accordingly SD_ri = 6*0.6161202 = 3.696721
  df$CESD_ri[1:10000] <- rpois(10000, 3.696721) # poisson as CESD mustn't be negative
  df$CESD_ri[10001:15500] <- df$CESD_ri[1:5500]
  # the SD for GAD7 is 3.41 as reported in Löwe et al.
  # accordingly SD_ri = 3.41*0.6161202 = 2.10097
  df$GAD_ri[1:10000] <- rpois(10000, 2.10097) # poisson as GAD7 mustn't be negative
  df$GAD_ri[10001:15500] <- df$GAD_ri[1:5500]
  # the SD for WMHV is 5384.3 as reported in Rodriguez et al.
  # accordingly SD_ri = 5384.3*0.6161202 = 3317.376 - this seems like a very disproportionate random intercept SD
  # to align the value more with the other variables we will half it 3317.376/2 = 1658.688
  df$WMHV_ri[1:10000] <- rpois(10000, 1658.688) # poisson as GAD7 mustn't be negative
  df$WMHV_ri[10001:15500] <- df$WMHV_ri[1:5500]
  
  # simulate outcome variables under the assumption of a linear relationship based on coefficients from Lammer et al.
  df$HCV <- as.vector(df$HCV_ri + 3828.409 + rnorm(15500,0,47.4) + df$age*(-5.697) + df$gender*(-46.390) + df$LSNS*(-5.697) +
                        df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
                        df$PSQI*(-11.915) + df$IPAQ*(-12.253))
  df$exfunct <- as.vector(df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
                            df$LSNS*(-0.017) +   df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
                            df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035))
  df$memo <- as.vector(df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
                         df$LSNS*(-0.008) +   df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
                         df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227))
  df$procspeed <- as.vector(df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
                              df$LSNS*(-0.017) +   df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
                              df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041))
  # to approximate the right-skewed distribution of CESD we will use a poission distribution 
  # to simulate the noise, the SD of the poisson distribution will be the residual SD calculated below 
  # to approximate the residuals' SD for these variables we will use the average residual SD of the four variables above relative to the variables SDs
  # HCV: residual sd = 47.4; ef: residual sd = 0.466906; memo: residual sd = 0.37214; procspeed: residual sd = 0.3964443
  # (0.466906 + 0.37214 + 0.3964443 + (47.4 / 420)) / 4 = 0.3370869
  # the CESD SD: 6 --> 6 * 0.3370869 = 2.022521
  # to get somewhat realistic values our scores should be similar to those from Lammer et al.
  # to get this: mean from Lammer et al. = 10 = mean (random intercept) + mean (noise) + intercept
  # 10 = 3.696721 + 2.022521 + 4.280758
  # Tang et al. regressed CESD-scores on LSNS - we will use the regression coefficients from their model 2a for LSNS, age and gender
  # they used a different version of the CESD, though - we will thus adjust the coefficients to the respective mean scores
  # CESD mean from Lammer et al. = 10
  # CESD mean from Tang et al. = 4.58
  # correction factor = 10 / 4.58 = 2.183406
  # coefficient for age from Tang et al.: 0.065 --> 0.065*2.183406 = 0.1419214
  # coefficient for gender from Tang et al.: 0.065 (negated as our default is female)  --> 0.065*2.183406 = 0.1419214
  # coefficient for LSNS from Tang et al.: 0.1145 (the average of the family and friends subscale, negated as our LSNS are inverted) 
  # --> 0.1145*2.183406 = 0.25
  # use absolute values for the rare case that neg values should be simulated
  df$CESD_out <- abs(round(df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.25, 0))
  
  # to approximate the right-skewed distribution of GAD7 we will use a poission distribution 
  # to simulate the noise: 0.3370869 * 3.41 = 1.149466 ---- 3.41 is SD from Löwe et al.
  # mean from Löwe et al. = 2.95 - already smaller than 1.149466 + 2.10097 --> no additional value added (subtracting could lead to neg. values and should thus be avoided)
  # Howren et al. regressed GAD7 scores on social isolation (based on LSNS)
  # they also regressed PHQ-9 scores (depressive symptoms) on social isolation (based on LSNS)
  # they used dichotomous isolated/not isolated predictor, though
  # to approximate a per LSNS point effect size we will use the effect size for anxiety relative to the effect size for depression from Howren et al.
  # the means of depression and anxiety scores are also different - thus we must also include this information in the calculation (we will use the baseline means)
  # and multiply it with the effect size from Tang et al. for depression
  # in the adjusted model the effect for anxiety is 1.91 (GAD7 mean = 7.93) and 2.24 for depression (PHQ-p mean = 9.36)
  # the relative effect size is: (1.91/7.93) / (2.24/9.36) = 1.00644
  # the absolute coefficient is: 0.25 * 1.00644 = 0.25161
  # Howren et al. controlled for gender and age but did not report their coefficients
  # as we are not too concerned about these variables for this simulation we will just use the coefficients for CESD adapted to the GAD7 mean
  # coefficient for age: 0.065 * (3.41/10) = 0.022165
  # coefficient for gender: 0.065 * (3.41/10) = 0.022165 # unlikely to be true but inconsequential for the simulation
  # modified coefficient for LSNS from  Howren et al.:  0.25161 
  # use absolute values for the rare case that neg values should be simulated
  df$GAD <- abs(round(df$GAD_ri + rpois(15500, 1.149466) + df$age*0.022165 + df$gender*0.022165 + df$LSNS*0.25161, 0))
  
  # to approximate the right-skewed distribution of WMHV we will use a poission distribution 
  # to simulate the noise: 0.3370869 * 5384.3 = 1814.977 ---- 5384.3 is SD from Rodriguez et al.
  # mean from Rodriguez et al. = 2952.6 ---- 2952.6 already smaller than 1814.977 + 1658.688 --> no additional value added (subtracting could lead to neg. values and should thus be avoided)
  # to approximate the coefficient we will additionlly draw information from Lammer et al.
  # in Rodriguez et al. the coefficient for social isolation (categorical) on log(WMHV) is 0.09
  # exp(0.09) = 1.094174
  # this can be interpreted as a 9% difference in WMHV in socially isolated individuals 
  # at the average WMHV of 2952.6 mm³ this would be an increase of 2952.6*0.094174 = 278.0582
  # this would be the effect size for categorical social isolation
  # to approximate the effect size per point on the LSNS we can use a sensitivity analysis from Lammer et al. (Appendix table 16) 
  # while the effect on HCV per point was -5.7, the effect of categorical LSNS was -81.625 --> -81.625 / -5.7 = 14.32018
  # accordingly our predictor for WMHV is 278.0582 / 14.32018 = 19.41723
  # we can approximately convert the predictors for age and gender into absolute values on a linear WMHV scale 
  # age: 0.05 --> exp(0.05) = 1.051271 --> 2952.6*0.051271 = 151.3828
  # gender: 0.30 --> exp(0.3) = 1.349859 --> 2952.6*0.349859 = 1032.994 (negate because our default is female)
  # use absolute values for the rare case that neg values should be simulated
  df$WMHV <- abs(round(df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 + df$LSNS*19.41723, 0))
  
  
  # CAVE: LSNS scores will be modified so that higher scores imply greater social isolation
  # we can simulate 4 scenarios in addition to the linear case
  # scenario 1: an additional exponential effect at LSNS > 18
  # scenario 2: an additional linear effect at LSNS > 18
  # scenario 3: an exponential effect only at LSNS > 18
  # scenario 4: a linear effect only at LSNS > 18
  # let's now simulate our outcome variables under the assumption of an additional exponential effect at LSNS > 18 (scenario 1)
  # the average slope across the scores 0-30 should be equivalent to the average slope under the linear assumption
  # for the threshold to have a substantial effect, we want the average slope above the threshold to be 3 times as large as below
  # for HCV:
  # linear: HCV = -5.697*LSNS + rest --> -5.697 = x
  # exponential for scores 0-18: HCV_exp = y*LSNS + rest
  # exponential for scores 19-30: HCV_exp = y*LSNS + z*((LSNS-18)^2) + rest
  # we can get the slope at each possible LSNS score by using the differential equations
  # in the linear case the slope is -5.697 at all LSNS-values
  # in the exponential case the slope is y for values < 19
  # in the exponential case the slope is y + 2z*(LSNS-18) for values > 18
  # to obtain the average slope, the slopes at each LSNS score should be weighted by their prevalence
  # the prevalence of each score can be derived from the density function of the normal distribution
  # we find that ca. 80% of observations are below threshold and 20% above threshold
  # sum(dnorm(min(df$LSNS):18, mean = 14.1, sd = 5.2)) = 0.8015948
  # sum(dnorm(19:max(df$LSNS), mean = 14.1, sd = 5.2)) = 0.1983501
  # thus, if we want to calculate the average slope in the exponential case across LSNS scores:
  # sum((1:(max(df$LSNS)-18))*dnorm(19:max(df$LSNS), mean = 14.1, sd = 5.2)) / sum(dnorm(19:max(df$LSNS), mean = 14.1, sd = 5.2))
  # = 3.407734
  # the average slope for LSNS > 18 = y + 3.4*2*z = y + 6.8z
  # x = 0.8*y + 0.2*(y+6.8z) = 0.8y + 0.2y + 0.2*6.8z = y + 1.36z
  # for the threshold to have a substantial effect, we want the average slope above the threshold to be 3 times as large as below
  # 3y = y+1.36z
  # 2y = 1.36z
  # y = 0.68z
  # x = -5.697 = y + 1.36z = 0.68z + 1.36z = 2.04z 
  # -5.697 / 2.04 = z = -2.792647
  # y = 0.68z = 0.68*(-2.792647) = -1.899
  
  # analogously the following coefficients for the other outcomes can be calculated
  # ef: x = -0.017;  z = -0.017 / 2.04 = -0.008333333; y = 0.68*(-0.008333333) = -0.005666666
  # memo: x = -0.008; z = -0.008 / 2.04 = -0.003921569; y = 0.68*(-0.003921569) = -0.002666667
  # procspeed: x = -0.017; z = -0.017 / 2.04 = -0.008333333; y = 0.68*(-0.008333333) = -0.005666666
  # CESD: x = 0.25; z = 0.25 / 2.04 = 0.122549; y = 0.68*0.122549 = 0.08333334
  # GAD7: x = 0.25161; z = 0.25161 / 2.04 = 0.1233382; y = 0.68*0.1233382 = 0.08387002
  # WMHV: x = 19.41723; z = 19.41723 / 2.04 = 9.51825; y = 0.68 * 9.51825 = 6.472411
  
  # simulate HCV for all subjects only using y
  df$HCVexp <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + df$LSNS*(-1.899) +
    df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
    df$PSQI*(-11.915) + df$IPAQ*(-12.253)
  # add exponential component if LSNS > 18
  df$HCVexp <- ifelse(df$LSNS > 18, df$HCVexp - 2.792647*((df$LSNS-18)^2), df$HCVexp)
  
  # do the same for other outcomes
  df$exfunctexp <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
    df$LSNS*(-0.005666666) +   df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
    df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
  df$exfunctexp <- ifelse(df$LSNS > 18, df$exfunctexp - 0.008333333*((df$LSNS-18)^2), df$exfunctexp)
  df$memoexp <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
    df$LSNS*(-0.002666667) +   df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
    df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
  df$memoexp <- ifelse(df$LSNS > 18, df$memoexp - 0.003921569*((df$LSNS-18)^2), df$memoexp)
  df$procspeedexp <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
    df$LSNS*(-0.005666666) +   df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
    df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)
  df$procspeedexp <- ifelse(df$LSNS > 18, df$procspeedexp - 0.008333333*((df$LSNS-18)^2), df$procspeedexp)
  df$CESD_outexp <- df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.08333334
  df$CESD_outexp <- ifelse(df$LSNS > 18, abs(round(df$CESD_outexp + 0.122549*((df$LSNS-18)^2),0)), abs(round(df$CESD_outexp,0)))
  df$GADexp <- df$GAD_ri + rpois(15500, 1.149466) + df$age*0.022165 + df$gender*0.022165 + df$LSNS*0.08387002
  df$GADexp <- ifelse(df$LSNS > 18, abs(round(df$GADexp + 0.1233382*((df$LSNS-18)^2),0)), abs(round(df$GADexp,0)))
  df$WMHVexp <- df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 + df$LSNS*6.472411
  df$WMHVexp <- ifelse(df$LSNS > 18, abs(round(df$WMHVexp + 9.51825*((df$LSNS-18)^2),0)), abs(round(df$WMHVexp,0)))
  
  # similarly we can simulate our outcome variables for the case that there is an additional linear effect of LSNS above threshold
  # the average slope across the scores 0-30 should be equivalent to the average slope under the linear assumption
  # we can get the slope at each possible LSNS score by using the differential equations
  # for HCV:
  # in the linear case the slope is -5.697 at all LSNS-values
  # in this scenario the slope is y for values < 19 and y + z for values > 18
  # -5.697 = y + 0.2*z
  # again the slope should be 3 times as large
  # 3y = y + z --> 2y = z
  # -5.697 = y + 0.2*2*y = 1.4y
  # -4.069286 = y
  # 2 * (-4.069286) = z = -8.138572
  df$HCVlin <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + 
    df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
    df$PSQI*(-11.915) + df$IPAQ*(-12.253) - 4.069286*df$LSNS
  # add linear component if LSNS > 18
  df$HCVlin <- ifelse(df$LSNS > 18, df$HCVlin - (2*4.069286)*(df$LSNS-18), df$HCVlin)
  
  # analogously the following coefficients for the other outcomes can be calculated
  # ef: x = -0.017;  y = -0.017 / 1.4 = -0.01214286; z = 2 * y = -0.02428572
  # memo: x = -0.008; y = -0.008 / 1.4 = -0.005714286; z = 2 * y = -0.01142857
  # procspeed: x = -0.017; y = -0.017 / 1.4 = -0.01214286; z = 2 * y = -0.02428572
  # CESD: x = 0.25; y = 0.25/1.4 = 0.1785714; z = 2 * y = 0.3571428
  # GAD7: x = 0.25161; y = 0.25161/1.4 = 0.1797214; z = 2 * y = 0.3594428
  # WMHV: x = 19.41723; y = 19.41723 / 1.4 = 13.86945; z = 2 * y = 27.7389
  
  df$exfunctlin <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
    df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
    df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035) - 0.01214286*df$LSNS
  df$exfunctlin <- ifelse(df$LSNS > 18, df$exfunctlin - (2*0.01214286)*(df$LSNS-18), df$exfunctlin)
  df$memolin <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
    df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
    df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227) - 0.005714286*df$LSNS
  df$memolin <- ifelse(df$LSNS > 18, df$memolin - (2*0.005714286)*(df$LSNS-18), df$memolin)
  df$procspeedlin <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
    df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
    df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041) - 0.01214286*df$LSNS
  df$procspeedlin <- ifelse(df$LSNS > 18, df$procspeedlin - (2*0.01214286)*(df$LSNS-18), df$procspeedlin)
  df$CESD_outlin <- df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.1785714
  df$CESD_outlin <- ifelse(df$LSNS > 18, abs(round(df$CESD_outlin + (2*0.1785714)*(df$LSNS-18),0)), abs(round(df$CESD_outlin,0)))
  df$GADlin <- df$GAD_ri + rpois(15500, 1.149466) + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.1797214
  df$GADlin <- ifelse(df$LSNS > 18, abs(round(df$GADlin + (2*0.1797214)*(df$LSNS-18),0)), abs(round(df$GADlin,0)))
  df$WMHVlin <- df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 + df$LSNS*13.86945
  df$WMHVlin <- ifelse(df$LSNS > 18, abs(round(df$WMHVlin + (2*13.86945)*(df$LSNS-18),0)), abs(round(df$WMHVlin,0)))
  
  # we shall now simulate our outcome variables for the thresholded exponential case (scenario 3)
  # for HCV: 
  # x = -5.697 
  # slope for LSNS < 19: 0
  # slope for LSNS > 18: 2z*(LSNS-18) 
  # -5.697 = 0.8*0 + 0.2*2*3.407734*z = 1.363094z
  # -4.179462 = z
  
  df$HCVthreshexp <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + 
    df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
    df$PSQI*(-11.915) + df$IPAQ*(-12.253)
  # add exponential component if LSNS > 18
  df$HCVthreshexp <- ifelse(df$LSNS > 18, df$HCVthreshexp - 4.179462*((df$LSNS-18)^2), df$HCVthreshexp)
  
  # analogously the following coefficients for the other outcomes can be calculated
  # ef: x = -0.017;  z = -0.017 / 1.363094 = -0.01247163
  # memo: x = -0.008; z = -0.008 / 1.363094 = -0.005869001
  # procspeed: x = -0.017; z = -0.017 / 1.363094 = -0.01247163
  # CESD: x = 0.25; z = 0.25 / 1.363094 = 0.1834063
  # GAD7: x = 0.25161; z = 0.25161 / 1.363094 = 0.1845874
  # WMHV: x = 19.41723 ; z = 19.41723  / 1.363094 = 14.24497 
  
  df$exfunctthreshexp <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
    df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
    df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
  df$exfunctthreshexp <- ifelse(df$LSNS > 18, df$exfunctthreshexp - 0.01247163*((df$LSNS-18)^2), df$exfunctthreshexp)
  df$memothreshexp <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
    df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
    df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
  df$memothreshexp <- ifelse(df$LSNS > 18, df$memothreshexp - 0.005869001*((df$LSNS-18)^2), df$memothreshexp)
  df$procspeedthreshexp <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
    df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
    df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)
  df$procspeedthreshexp <- ifelse(df$LSNS > 18, df$procspeedthreshexp - 0.01247163*((df$LSNS-18)^2), df$procspeedthreshexp)
  df$CESD_outthreshexp <- df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214
  df$CESD_outthreshexp <- ifelse(df$LSNS > 18, abs(round(df$CESD_outthreshexp + 0.1834063*((df$LSNS-18)^2),0)), abs(round(df$CESD_outthreshexp,0)))
  df$GADthreshexp <- df$GAD_ri + rpois(15500, 1.149466) + df$age*0.1419214 + df$gender*0.1419214
  df$GADthreshexp <- ifelse(df$LSNS > 18, abs(round(df$GADthreshexp + 0.1845874*((df$LSNS-18)^2),0)), abs(round(df$GADthreshexp,0)))
  df$WMHVthreshexp <- df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 
  df$WMHVthreshexp <- ifelse(df$LSNS > 18, abs(round(df$WMHVthreshexp + 14.24497*((df$LSNS-18)^2),0)), abs(round(df$WMHVthreshexp,0)))
  
  # lastly, we simulate our outcome variables for the thresholded linear case (scenario 4)
  # for HCV: 
  # x = -5.697 
  # slope for LSNS < 19: 0
  # slope for LSNS > 18: y 
  # -5.697 = 0.8+0 + 0.2*y = 0.2y
  # -28.485 = y
  
  df$HCVthreshlin <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + 
    df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
    df$PSQI*(-11.915) + df$IPAQ*(-12.253)
  # add linear component if LSNS > 18
  df$HCVthreshlin <- ifelse(df$LSNS > 18, df$HCVthreshexp - 28.485*(df$LSNS-18), df$HCVthreshlin)
  
  # analogously for other outcomes
  # ef: x = -0.017;  y = -0.017*5 = -0.085
  # memo: x = -0.008; y = -0.008*5 = -0.04
  # procspeed: x = -0.017; y = -0.017*5 = -0.085
  # CESD: x = 0.25; y = 0.25*5 = 1.25
  # GAD7: x = 0.25161; y = 0.25161*5 = 1.25805
  # WMHV: x = 19.41723; y = 19.41723 * 5  = 97.08615
  
  df$exfunctthreshlin <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
    df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
    df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
  df$exfunctthreshlin <- ifelse(df$LSNS > 18, df$exfunctthreshlin - 0.085*(df$LSNS-18), df$exfunctthreshlin)
  df$memothreshlin <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
    df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
    df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
  df$memothreshlin <- ifelse(df$LSNS > 18, df$memothreshlin - 0.04*(df$LSNS-18), df$memothreshlin)
  df$procspeedthreshlin <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
    df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
    df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)
  df$procspeedthreshlin <- ifelse(df$LSNS > 18, df$procspeedthreshlin - 0.085*(df$LSNS-18), df$procspeedthreshlin)
  df$CESD_outthreshlin <- df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214
  df$CESD_outthreshlin <- ifelse(df$LSNS > 18, abs(round(df$CESD_outthreshlin + 1.25*(df$LSNS-18),0)), abs(round(df$CESD_outthreshlin,0)))
  df$GADthreshlin <- df$GAD_ri + rpois(15500, 1.149466) + df$age*0.1419214 + df$gender*0.1419214
  df$GADthreshlin <- ifelse(df$LSNS > 18, abs(round(df$GADthreshlin + 1.25805*(df$LSNS-18),0)), abs(round(df$GADthreshlin,0)))
  df$WMHVthreshlin <- df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 
  df$WMHVthreshlin <- ifelse(df$LSNS > 18, abs(round(df$WMHVthreshlin + 97.08615*(df$LSNS-18),0)), abs(round(df$WMHVthreshlin,0)))
  
  # run models
  for (var in outvars_gauss){
    simres[simres$simulation == n & simres$outcome == var,c("AICbase", "BICbase", "pbase", "AICfull", "BICfull", "pfull", "edf")] <- simulate_amm(dv = var)
  }
  for (var in outvars_poisson){
    simres[simres$simulation == n & simres$outcome == var,c("AICbase", "BICbase", "pbase", "AICfull", "BICfull", "pfull", "edf")] <- simulate_gamm(dv = var)
  }
  save.image("Simulation.RData")
}

simres$dAICbasefull <- simres$AICbase - simres$AICfull
simres$dBICbasefull <- simres$BICbase - simres$BICfull
simres_summary <- data.frame(matrix(ncol = 12, nrow = 5*7))
colnames(simres_summary) <- c("outcome", "meanAICbase", "meanBICbase", "meanAICfull", "meanBICfull", "meanpfull", "propsignif", "meandAICbasefull", 
                              "propAICprofull", "meandBICbasefull", "propBICprofull", "meanedf")
simres_summary$outcome <- as.vector(outer(c("HCV", "exfunct", "memo", "procspeed", "CESD_out", "WMHV", "GAD"),c("", "exp", "lin", "threshexp", "threshlin"), paste, sep = "" ))
for (outcome in simres_summary$outcome){
  for (column_name in c("AICbase", "BICbase", "AICfull", "BICfull", "pfull", "dAICbasefull", "dBICbasefull", "edf")) {
    simres_summary[simres_summary$outcome == outcome,c(paste0("mean", column_name))] <- mean(simres[simres$outcome == outcome,column_name], na.rm = T)
  }
  simres_summary[simres_summary$outcome == outcome,"propAICprofull"] <- length(which(simres[simres$outcome == outcome,"dAICbasefull"] > 0)) / length(which(!is.na(simres[simres$outcome == outcome,"dAICbasefull"])))
  simres_summary[simres_summary$outcome == outcome,"propBICprofull"] <- length(which(simres[simres$outcome == outcome,"dBICbasefull"] > 0)) / length(which(!is.na(simres[simres$outcome == outcome,"dBICbasefull"])))
  simres_summary[simres_summary$outcome == outcome,"propsignif"] <- length(which(simres[simres$outcome == outcome,"pfull"] < 0.05)) / length(which(!is.na(simres[simres$outcome == outcome,"pfull"])))
}
write.csv(simres_summary, file = "simresults.csv", quote = F, row.names = F)
save.image("Simulation.RData")

