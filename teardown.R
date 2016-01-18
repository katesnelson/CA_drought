
library(memisc)
library(dplyr)
library (lme4)
library (lmerTest)


#setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/model1")

dsf <- tbl_df(read.csv("./pt_data_kn.txt", stringsAsFactors = FALSE))
ds = ds[ds$farm_flag==1,]

#change datatype 
ds = transform(ds, 
          tvp_09 = as.numeric(tvp_09),
          tvp_14 = as.numeric(tvp_14),
          delta_tvp = as.numeric(delta_tvp),
          AreaAcres = as.numeric(AreaAcres))

#change values of 999.99 to missing value (Perc_Rip and Perc_Pre1914)
ds[ds==999.99] <- NA

#null model using FIML estimation
null = lmer(delta_tvp ~ 1 + (1 | HUC12), data = ds, REML = "False")
summary (null)

#calculate variance components and extract group variance and residual variance
vc = VarCorr(null)
tau = as.data.frame(vc)[1,4]
resid = as.data.frame(vc)[2,4]

# calculate ICC
totvar = tau + resid
ICC = tau /totvar
ICC

#full model
GW_density_change = ds$GW_dnsty15 - ds$GW_dnsty09
lc.model = lmer(delta_tvp ~ farm_flag + delta_lc + diverse + WR_density + GW_dnsty15 + GW_dnsty09 + (1 | HUC12), data = ds, REML = "False") 
summary(lc.model)

#start with all effects as random; if significant keep, if not drop
#remove RE from slopes before intercepts
#drop L2 before L1
#any cross level interactions?
#compare AIC at each step for improvements in model fit






devcomp = getME(lc.model,"devcomp")
dev.lc = as.numeric(devcomp$cmp[8])
dev.lc
devstat = dev -dev.lc
devstat #not significant, regardless of what q is

#add landuse change as predictor of int and slope
lc2.model = lmer(delta_tvp ~ delta_lc + (delta_lc | HUC12), data = ds, REML = "False")
summary(lc2.model)
devcomp = getME(lc2.model,"devcomp")
dev.lc2 = as.numeric(devcomp$cmp[8])
dev.lc2
devstat = dev -dev.lc2
devstat # difference in number of est parameters (q) should be 2, this is significant
