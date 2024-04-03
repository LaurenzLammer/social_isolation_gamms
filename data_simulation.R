library(lmerTest)

simres <- data.frame(matrix(nrow = 500, ncol = 6))
colnames(simres) <- c("simulation", "outcome", "coeff thresh", "p thresh", "coeff exp", "p exp")
simres$simulation <- unlist(lapply(1:100, rep, times = 5))
simres$outcome <- rep(c("HCV", "HCVexp", "HCVlin", "HCVthreshexp", "HCVthreshlin"))

df <- data.frame(matrix(ncol = 41, nrow = 10000+5500))
colnames(df) <- c("ID", "fu", "neurocog", "LSNS", "CESD", "GAD", "WMHV", "HCV", "memo", "procspeed", 
                  "exfunct", "CESDexp", "GADexp", "WMHVexp", "HCVexp", "memoexp", "procspeedexp", 
                  "exfunctexp", "CESDthresh", "GADthresh", "WMHVthresh", "HCVthresh", "memothresh", "procspeedthresh", 
                  "exfunctthresh", "HCV_ri", "WMHV_ri", "memo_ri", "procspeed_ri", "exfunct_ri", "CESD_ri", "GAD_ri", "BMI", 
                  "hypertension", "diabetes", "education", "age", "gender", "SES", "PSQI", "IPAQ")
df$ID[1:10000] <- as.character(1:10000)
df$ID[10001:15500] <- as.character(1:5500)
df$fu[1:10000] <- 0
df$fu[10001:15500] <- 1
df$neurocog[c(1:2500, 10001:11800)] <- 1
df$neurocog[-c(1:2500, 10001:11800)] <- 0

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
df$PSQI <- scale(df$PSQI)
df$IPAQ <- scale(df$IPAQ)
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
df$HCV_ri[1:10000] <- rnorm(10000, 0, 309.07)
df$HCV_ri[10001:15500] <- df$HCV_ri[1:5500]
#df$WMHV_ri[1:10000] <- rnorm(10000, 0, XXX)
#df$WMHV_ri[10001:15500] <- df$WMHV_ri[1:5500]
df$exfunct_ri[1:10000] <- rnorm(10000, 0, 0.6415)
df$exfunct_ri[10001:15500] <- df$exfunct_ri[1:5500]
df$memo_ri[1:10000] <- rnorm(10000, 0, 0.5878)
df$memo_ri[10001:15500] <- df$memo_ri[1:5500]
df$procspeed_ri[1:10000] <- rnorm(10000, 0, 0.4993)
df$procspeed_ri[10001:15500] <- df$procspeed_ri[1:5500]
#df$CESD_ri[1:10000] <- rnorm(10000, 0, 309.07)
#df$CESD_ri[10001:15500] <- df$CESD_ri[1:5500]
#df$GAD_ri[1:10000] <- rnorm(10000, 0, 309.07)
#df$GAD_ri[10001:15500] <- df$GAD_ri[1:5500]

df$HCV <- df$HCV_ri + 3828.409 + rnorm(15500,0,47.4) + df$age*(-5.697) + df$gender*(-46.390) + df$LSNS*(-5.697) +
  df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
  df$PSQI*(-11.915) + df$IPAQ*(-12.253)
df$exfunct <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
  df$LSNS*(-0.017) +   df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
  df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
df$memo <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
  df$LSNS*(-0.008) +   df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
  df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
df$procspeed <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
  df$LSNS*(-0.017) +   df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
  df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)

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
# exponential for scores 19-30: HCV_exp = y + z*((LSNS-18)^2) + rest
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

# analogously the following coefficients for the cognitive functions can be calculated
# ef: x = -0.017;  z = -0.017 / 10.2232 = -0.001662884; y = 3.407734*(-0.001662884) = -0.005666666
# memo: x = -0.008; z = -0.008 / 10.2232 = -0.0007825338; y = 3.407734*(-0.0007825338) = -0.002666667
# procspeed: x = -0.017; z = -0.017 / 10.2232 = -0.001662884; y = 3.407734*(-0.001662884) = -0.005666666


# simulate HCV for all subjects only using y
df$HCVexp <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + df$LSNS*(-1.899) +
  df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
  df$PSQI*(-11.915) + df$IPAQ*(-12.253)
# add exponential component if LSNS > 18
df$HCVexp <- ifelse(df$LSNS > 18, df$HCVexp - 0.5572619*((df$LSNS-18)^2), df$HCVexp)

# do the same for cog functs
df$exfunctexp <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
  df$LSNS*(-0.005666666) +   df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
  df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
df$exfunctexp <- ifelse(df$LSNS > 18, df$exfunctexp - 0.001662884*(df$LSNS-18), df$exfunctexp)
df$memoexp <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
  df$LSNS*(-0.002666667) +   df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
  df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
df$memoexp <- ifelse(df$LSNS > 18, df$memoexp - 0.0007825338*(df$LSNS-18), df$memoexp)
df$procspeedexp <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
  df$LSNS*(-0.005666666) +   df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
  df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)
df$procspeedexp <- ifelse(df$LSNS > 18, df$procspeedexp - 0.001662884*(df$LSNS-18), df$procspeedexp)

# similarly we can simulate our outcome variables for the case that there is an additional linear effect of LSNS above threshold
# the average slope across the scores 0-30 should be equivalent to the average slope under the linear assumption
# we can get the slope at each possible LSNS score by using the differential equations
# for HCV:
# in the linear case the slope is -5.697 at all LSNS-values
# in this scenario the slope is y for values < 19 and z for values > 18
# -5.697 = 0.8*y + 0.2*z
# again the slope should be 3 times as large
# 3y = z
# -5.697 = 0.8*y + 0.2*3*y = 1.4y
# -4.069286 = y
# 3 * (-4.069286) = z = -12.20786
df$HCVlin <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + 
  df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
  df$PSQI*(-11.915) + df$IPAQ*(-12.253) - 4.069286*df$LSNS
# add linear component if LSNS > 18
df$HCVlin <- ifelse(df$LSNS > 18, df$HCVlin - (12.20786 - 4.069286)*(df$LSNS-18), df$HCVlin)

# analogously the following coefficients for the cognitive functions can be calculated
# ef: x = -0.017;  y = -0.017 / 1.4 = -0.01214286; z = 3 * z = -0.03642858
# memo: x = -0.008; y = -0.008 / 1.4 = -0.005714286; y = 3 * z = -0.01714286
# procspeed: x = -0.017; y = -0.017 / 1.4 = -0.01214286; y = 3 * z = -0.03642858
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

# we shall now simulate our outcome variables for the thresholded exponential case (scenario 3)
# for HCV: 
# x = -5.697 
# slope for LSNS < 19: 0
# slope for LSNS > 18: 2z*(LSNS-18) 
# -5.697 = 0.8*0 + 0.2*3.407734*z = 0.6815468z
# -8.358927 = z

df$HCVthreshexp <- df$HCV_ri + 3828.409 + rnorm(n = 15500,mean = 0, sd = 47.4) + df$age*(-5.697) + df$gender*(-46.390) + 
  df$BMI*13.722 + df$CESD*12.184 + df$diabetes*(-99.613) + df$education*(-93.450) + df$hypertension*(-21.182) + 
  df$PSQI*(-11.915) + df$IPAQ*(-12.253)
# add exponential component if LSNS > 18
df$HCVthreshexp <- ifelse(df$LSNS > 18, df$HCVthreshexp - 8.358927*((df$LSNS-18)^2), df$HCVthreshexp)

# analogously the following coefficients for the cognitive functions can be calculated
# ef: x = -0.017;  z = -0.017 / 0.6815468 = -0.02494326
# memo: x = -0.008; z = -0.008 / 0.6815468 = -0.01173801
# procspeed: x = -0.017; z = -0.017 / 0.6815468 = -0.02494326

df$exfunctthreshexp <- df$exfunct_ri + 0.643381 + rnorm(15500,0,0.466906) + df$age*(-0.014) + df$gender*(-0.122) + 
  df$BMI*(-0.074) + df$CESD*(-0.141) + df$diabetes*(-0.055) + df$education*(-0.331) + 
  df$hypertension*(-0.097) + df$PSQI*(0.024) + df$IPAQ*(-0.035)
df$exfunctthreshexp <- ifelse(df$LSNS > 18, df$exfunctthreshexp - 0.02494326*((df$LSNS-18)^2), df$exfunctthreshexp)
df$memothreshexp <- df$memo_ri + 0.479806 + rnorm(15500,0,0.37214) + df$age*(-0.033) + df$gender*(-0.421) + 
  df$BMI*(-0.031) + df$CESD*(-0.120) + df$diabetes*(-0.038) + df$education*(-0.173) + 
  df$hypertension*(0.017) + df$PSQI*(-0.0046) + df$IPAQ*(-0.0227)
df$memothreshexp <- ifelse(df$LSNS > 18, df$memothreshexp - 0.01173801*((df$LSNS-18)^2), df$memothreshexp)
df$procspeedthreshexp <- df$procspeed_ri + 0.462252 + rnorm(15500,0,0.3964443) + df$age*(-0.036) + df$gender*(-0.124) + 
  df$BMI*(-0.013) + df$CESD*(-0.027) + df$diabetes*(-0.005) + df$education*(-0.115) + 
  df$hypertension*(-0.067) + df$PSQI*(-0.0566) + df$IPAQ*(-0.041)
df$procspeedthreshexp <- ifelse(df$LSNS > 18, df$procspeedthreshexp - 0.02494326*((df$LSNS-18)^2), df$procspeedthreshexp)

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

# analogously for cog functs
# ef: x = -0.017;  y = -0.017*5 = -0.085
# memo: x = -0.008; y = -0.008*5 = -0.04
# procspeed: x = -0.017; y = -0.017*5 = -0.085

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

