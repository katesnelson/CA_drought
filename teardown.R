
library(memisc)
library(dplyr)
library (lme4)
library (lmerTest)

dsf <- tbl_df(read.csv("C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Kate_Emily\\Data\\pt_data_kn.txt", stringsAsFactors = FALSE))
ds = dsf[dsf$farm_flag==1,]

#change datatype 
ds = transform(ds, 
          tvp_09 = as.numeric(tvp_09),
          tvp_14 = as.numeric(tvp_14),
          delta_tvp = as.numeric(delta_tvp),
          AreaAcres = as.numeric(AreaAcres))

#change values of 999.99 to missing value (Perc_Rip and Perc_Pre1914)
ds[ds==999.99] <- NA

#farming pixels only
data <- ds[ds$farm_flag == 1,]

#null model using FIML estimation
null = lmer(delta_tvp ~ 1 + (1 | HUC12), data = data, REML = "False")
summary (null)

#calculate variance components and extract group variance and residual variance
vc = VarCorr(null)
tau = as.data.frame(vc)[1,4]
resid = as.data.frame(vc)[2,4]
# calculate ICC
totvar = tau + resid
ICC = tau /totvar
ICC

#non-mlm
GLM <- glm(delta_tvp ~ delta_lc + diverse + Avg_WSEL_5yrChange + WR_density + Perc_Rip + Perc_Pre1914 + GW_dnsty15, data = ds)
summary(GLM)


#http://www.rensenieuwenhuis.nl/r-sessions-16-multilevel-model-specification-lme4/
#full model, varying intercepts and slopes for all variables
#no reason to think the effects of land use change should vary across HUCs
#reason to think that the effects of diverse, WR_density and GW_dnsty15 could vary
M1 = lmer(delta_tvp ~ delta_lc + diverse + Avg_WSEL_5yrChange + WR_density + Perc_Rip + Perc_Pre1914 + GW_dnsty15 + 
          delta_lc*diverse + delta_lc*WR_density +delta_lc*GW_dnsty15 + (1 + delta_lc|HUC12), data = ds, REML = "False") 

summary(M1)
devcomp = getME(M1,"devcomp")
devM1 = as.numeric(devcomp$cmp[8])

#drop perc_pre1914
M2 = lmer(delta_tvp ~ delta_lc + diverse + Avg_WSEL_5yrChange + WR_density + Perc_Rip + GW_dnsty15 + 
            delta_lc*diverse + delta_lc*WR_density +delta_lc*GW_dnsty15 + (1 + delta_lc|HUC12), data = ds, REML = "False") 

summary(M2)
devcomp = getME(M2,"devcomp")
devM2 = as.numeric(devcomp$cmp[8])
devstat = devM2 - devM1
devstat 

#drop Avg_WSEL_5yrChange
M3 = lmer(delta_tvp ~ delta_lc + diverse + WR_density + Perc_Rip + GW_dnsty15 + Perc_Pre1914 +
            delta_lc*diverse + delta_lc*WR_density +delta_lc*GW_dnsty15 + (1 + delta_lc|HUC12), data = ds, REML = "False") 

summary(M3)
devcomp = getME(M3,"devcomp")
devM3 = as.numeric(devcomp$cmp[8])
devstat = devM3 - devM1
devstat 

#drop diverse*delat_LC
M4 = lmer(delta_tvp ~ delta_lc + diverse + WR_density + Perc_Rip + GW_dnsty15 + Perc_Pre1914 + Avg_WSEL_5yrChange +
            delta_lc*WR_density +delta_lc*GW_dnsty15 + (1 + delta_lc|HUC12), data = ds, REML = "False") 

summary(M4)
devcomp = getME(M4,"devcomp")
devM4 = as.numeric(devcomp$cmp[8])
devstat = devM4 - devM1
devstat 

#drop diverse*delat_LC, delta_lc*WR_density
M5 = lmer(delta_tvp ~ delta_lc + diverse + WR_density + Perc_Rip + GW_dnsty15 + Perc_Pre1914 + Avg_WSEL_5yrChange +
            delta_lc*GW_dnsty15 + (1 + delta_lc|HUC12), data = ds, REML = "False") 

summary(M5)
devcomp = getME(M5,"devcomp")
devM5 = as.numeric(devcomp$cmp[8])
devstat = devM5 - devM4
devstat 

#drop diverse*delat_LC, delta_lc*WR_density, delta_lc*GW_dnsty15
M6 = lmer(delta_tvp ~ delta_lc + diverse + WR_density + Perc_Rip + GW_dnsty15 + Perc_Pre1914 + Avg_WSEL_5yrChange +
            (1 + delta_lc|HUC12), data = ds, REML = "False") 

summary(M6)
devcomp = getME(M6,"devcomp")
devM6 = as.numeric(devcomp$cmp[8])
devstat = devM6 - devM5
devstat 

#devstat_complex - dec_lesscomplex = +, that's a good sign
#depending on number of parameters that have changed, q = 1 and devstat < 3.8 then the improvement in fit is not significant 

