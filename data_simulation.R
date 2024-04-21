# load required package
library(lmerTest)
# define functions for LME modelling
simulate_lmer <- function(dv){ # dv is respective dependent variable
  mod <- lmer(data = df[df$neurocog == 1,], formula = paste0(dv, " ~ age + gender + LSNS + thresh + (1|ID)"))
  smod <- summary(mod)
  # save results into external dataframe
  return(smod$coefficients[c("thresh", "LSNS"), c("Estimate", "Pr(>|t|)")])
}
simulate_glmer <- function(dv){
  mod <- glmer(data = df, formula = paste0(dv, " ~ age + gender + LSNS + thresh + (1|ID)"), family = poisson(link = "log"))
  smod <- summary(mod)
  return(smod$coefficients[c("thresh", "LSNS"), c("Estimate", "Pr(>|z|)")])
}

# prepare data frame for simulation results
simres <- data.frame(matrix(nrow = 100*5*7, ncol = 8))
colnames(simres) <- c("simulation", "outcome", "coeffthresh", "pthresh", "pthresh_sided", "coeffLSNS", "p_LSNS", "p_LSNS_sided")
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
  # while the effect on HCV per point was –5.7, the effect of categorical LSNS was –81.625 --> -81.625 / -5.7 = 14.32018
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
  # x = 0.8*y + 0.2*(y+6.8z) = 0.8y + 0.2y + 6.8z = y + 6.8z
  # for the threshold to have a substantial effect, we want the average slope above the threshold to be 3 times as large as below
  # 3y = y+6.8z
  # 2y = 6.8z
  # y = 3.407734z
  # x = -5.697 = y + 6.8z = 3.407734z + 6.8z = 3 * 3.407734z = 10.2232z
  # -5.697 / 10.2232 = z = -0.5572619
  # y = 3.407734z = 3.407734*(-0.5572619) = -1.899
  
  # analogously the following coefficients for the other outcomes can be calculated
  # ef: x = -0.017;  z = -0.017 / 10.2232 = -0.001662884; y = 3.407734*(-0.001662884) = -0.005666666
  # memo: x = -0.008; z = -0.008 / 10.2232 = -0.0007825338; y = 3.407734*(-0.0007825338) = -0.002666667
  # procspeed: x = -0.017; z = -0.017 / 10.2232 = -0.001662884; y = 3.407734*(-0.001662884) = -0.005666666
  # CESD: x = 0.25; z = 0.25 / 10.2232 = 0.02445418; y = 3.407734*0.02445418 = 0.08333334
  # GAD7: x = 0.25161; z = 0.25161 / 10.2232 = 0.02461167; y = 3.407734*0.02461167 = 0.08387002
  # WMHV: x = 19.41723; z = 19.41723 / 10.2232 = 1.89933; y = 3.407734 * 1.89933 = 6.472411
  
  # simulate HCV for all subjects only using y
  df$HCVexp <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + df$LSNS*(-1.899) +
    df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
    df$PSQI*(-11.915) + df$IPAQ*(-12.253)
  # add exponential component if LSNS > 18
  df$HCVexp <- ifelse(df$LSNS > 18, df$HCVexp - 0.5572619*((df$LSNS-18)^2), df$HCVexp)
  
  # do the same for other outcomes
  df$exfunctexp <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
    df$LSNS*(-0.005666666) +   df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
    df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
  df$exfunctexp <- ifelse(df$LSNS > 18, df$exfunctexp - 0.001662884*((df$LSNS-18)^2), df$exfunctexp)
  df$memoexp <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
    df$LSNS*(-0.002666667) +   df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
    df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
  df$memoexp <- ifelse(df$LSNS > 18, df$memoexp - 0.0007825338*((df$LSNS-18)^2), df$memoexp)
  df$procspeedexp <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
    df$LSNS*(-0.005666666) +   df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
    df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)
  df$procspeedexp <- ifelse(df$LSNS > 18, df$procspeedexp - 0.001662884*((df$LSNS-18)^2), df$procspeedexp)
  df$CESD_outexp <- df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.08333334
  df$CESD_outexp <- ifelse(df$LSNS > 18, abs(round(df$CESD_outexp + 0.02445418*((df$LSNS-18)^2),0)), abs(round(df$CESD_outexp,0)))
  df$GADexp <- df$GAD_ri + rpois(15500, 1.149466) + df$age*0.022165 + df$gender*0.022165 + df$LSNS*0.08387002
  df$GADexp <- ifelse(df$LSNS > 18, abs(round(df$GADexp + 0.02461167*((df$LSNS-18)^2),0)), abs(round(df$GADexp,0)))
  df$WMHVexp <- df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 + df$LSNS*6.472411
  df$WMHVexp <- ifelse(df$LSNS > 18, abs(round(df$WMHVexp + 1.89933*((df$LSNS-18)^2),0)), abs(round(df$WMHVexp,0)))
  
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
  df$HCVlin <- ifelse(df$LSNS > 18, df$HCVlin - (12.20786 - 4.069286)*(df$LSNS-18), df$HCVlin)
  
  # analogously the following coefficients for the other outcomes can be calculated
  # ef: x = -0.017;  y = -0.017 / 1.4 = -0.01214286; z = 3 * y = -0.03642858
  # memo: x = -0.008; y = -0.008 / 1.4 = -0.005714286; z = 3 * y = -0.01714286
  # procspeed: x = -0.017; y = -0.017 / 1.4 = -0.01214286; z = 3 * y = -0.03642858
  # CESD: x = 0.25; y = 0.25/1.4 = 0.1785714; z = 3 * y = 0.5357142
  # GAD7: x = 0.25161; y = 0.25161/1.4 = 0.1797214; z = 3 * y = 0.5391642
  # WMHV: x = 19.41723; y = 19.41723 / 1.4 = 13.86945; z = 3 * y = 41.60835
  
  df$exfunctlin <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
    df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
    df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035) - 0.01214286*df$LSNS
  df$exfunctlin <- ifelse(df$LSNS > 18, df$exfunctlin - (0.036428586 - 0.01214286)*(df$LSNS-18), df$exfunctlin)
  df$memolin <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
    df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
    df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227) - 0.005714286*df$LSNS
  df$memolin <- ifelse(df$LSNS > 18, df$memolin - (0.01714286 -  - 0.005714286)*(df$LSNS-18), df$memolin)
  df$procspeedlin <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
    df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
    df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041) - 0.01214286*df$LSNS
  df$procspeedlin <- ifelse(df$LSNS > 18, df$procspeedlin - (0.036428586 - 0.01214286)*(df$LSNS-18), df$procspeedlin)
  df$CESD_outlin <- df$CESD_ri + rpois(15500, 2.022521) + 4.280758 + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.1785714
  df$CESD_outlin <- ifelse(df$LSNS > 18, abs(round(df$CESD_outlin + (0.5357142-0.1785714)*(df$LSNS-18),0)), abs(round(df$CESD_outlin,0)))
  df$GADlin <- df$GAD_ri + rpois(15500, 1.149466) + df$age*0.1419214 + df$gender*0.1419214 + df$LSNS*0.1797214
  df$GADlin <- ifelse(df$LSNS > 18, abs(round(df$GADlin + (0.5391642-0.1797214)*(df$LSNS-18),0)), abs(round(df$GADlin,0)))
  df$WMHVlin <- df$WMHV_ri + rpois(15500, 1814.977) + df$age*151.3828 + df$gender*1032.994 + df$LSNS*13.86945
  df$WMHVlin <- ifelse(df$LSNS > 18, abs(round(df$WMHVlin + (41.60835-13.86945)*(df$LSNS-18),0)), abs(round(df$WMHVlin,0)))
  
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
    simres[simres$simulation == n & simres$outcome == var,c("coeffthresh", "coeffLSNS", "pthresh", "p_LSNS")] <- simulate_lmer(dv = var)
  }
  for (var in outvars_poisson){
    simres[simres$simulation == n & simres$outcome == var,c("coeffthresh", "coeffLSNS", "pthresh", "p_LSNS")] <- simulate_glmer(dv = var)
  }
  
}
simres[simres$outcome %in% outvars_gauss,"pthresh_sided"] <- ifelse(simres[simres$outcome %in% outvars_gauss,"coeffthresh"] 
  < 0, simres[simres$outcome %in% outvars_gauss,"pthresh"]/2, 1-(simres[simres$outcome %in% outvars_gauss,"pthresh"]/2))
simres[simres$outcome %in% outvars_gauss,"p_LSNS_sided"] <- ifelse(simres[simres$outcome %in% outvars_gauss,"coeffLSNS"] 
  < 0, simres[simres$outcome %in% outvars_gauss,"p_LSNS"]/2, 1-(simres[simres$outcome %in% outvars_gauss,"p_LSNS"]/2))
simres[simres$outcome %in% outvars_poisson,"pthresh_sided"] <- ifelse(simres[simres$outcome %in% outvars_poisson,"coeffthresh"] 
  > 0, simres[simres$outcome %in% outvars_poisson,"pthresh"]/2, 1-(simres[simres$outcome %in% outvars_poisson,"pthresh"]/2))
simres[simres$outcome %in% outvars_poisson,"p_LSNS_sided"] <- ifelse(simres[simres$outcome %in% outvars_poisson,"coeffLSNS"] 
  > 0, simres[simres$outcome %in% outvars_poisson,"p_LSNS"]/2, 1-(simres[simres$outcome %in% outvars_poisson,"p_LSNS"]/2))

for (n in nrow(simres)){
  if (simres[n,"outcome"] %in% outvars_gauss){
    if (simres[n,"coeffthresh"] < 0){
      simres[n,"pthresh_sided"] <- simres[n,"pthresh"]*2
    } else {
      simres[n,"pthresh_sided"] <- 1-(simres[n,"pthresh"]/2)
    }
    if (simres[n,"coeffLSNS"] < 0){
      simres[n,"p_LSNS_sided"] <- simres[n,"p_LSNS"]*2
    } else {
      simres[n,"p_LSNS_sided"] <- 1-(simres[n,"p_LSNS"]/2)
    }
  }
  else {
    if (simres[n,"coeffthresh"] > 0){
      simres[n,"pthresh_sided"] <- simres[n,"pthresh"]*2
    } else {
      simres[n,"pthresh_sided"] <- 1-(simres[n,"pthresh"]/2)
    }
    if (simres[n,"coeffLSNS"] > 0){
      simres[n,"p_LSNS_sided"] <- simres[n,"p_LSNS"]*2
    } else {
      simres[n,"p_LSNS_sided"] <- 1-(simres[n,"p_LSNS"]/2)
    }
  }
}

simres_summary <- data.frame(matrix(ncol = 11, nrow = 5*7))
column_names <- c("outcome", "coeffthresh", "pthresh", "pthreshsign", "pthresh_sided", "pthresh_sidedsign", "coeffLSNS", 
                  "p_LSNS", "p_LSNSsign", "p_LSNS_sided", "p_LSNS_sidedsign")
colnames(simres_summary) <- column_names
simres_summary$outcome <- as.vector(outer(c("HCV", "exfunct", "memo", "procspeed", "CESD_out", "WMHV", "GAD"),c("", "exp", "lin", "threshexp", "threshlin"), paste, sep = "" ))
for (outcome in outcomes){
  for (column_name in c("coeffthresh", "pthresh", "pthresh_sided", "coeffLSNS", "p_LSNS", "p_LSNS_sided")) {
    simres_summary[simres_summary$outcome == outcome,column_name] <- mean(simres[simres$outcome == outcome,column_name])
  }
  for (column_name in c("pthreshsign", "pthresh_sidedsign", "p_LSNSsign", "p_LSNS_sidedsign")) {
    simres_summary[simres_summary$outcome == outcome,column_name] <- length(which(simres[simres$outcome == outcome,substr(column_name,1,nchar(column_name)-4)] < 0.05))
  }
}
false_positive_rate <- mean(simres_summary$pthresh_sidedsign[1:7]) # 6.142857
true_positive_rate <- mean(simres_summary$pthresh_sidedsign[8:35]) # 80.85714

