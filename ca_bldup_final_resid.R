
library(memisc)
library(dplyr)
library (lme4)
library (lmerTest)
library(HLMdiag)
library (gstat)
library(sp)
library(nlme)
library(spdep)

setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/model1")

ds <- tbl_df(read.csv("./pt_data_kn.txt", stringsAsFactors = FALSE))

#check data type
sapply(ds, class)

#change datatype 
ds = transform(ds, 
          tvp_09 = as.numeric(tvp_09),
          tvp_14 = as.numeric(tvp_14),
          delta_tvp = as.numeric(delta_tvp),
	    lc_09 = as.factor(lc_09),
          lc_14 = as.factor(lc_14),
          delta_lc = as.factor (delta_lc),
          AreaAcres = as.numeric(AreaAcres))
ds = mutate(ds, Rip = COUNT_Rip/fmmp_area, P1914 = COUNT_Pre1914/fmmp_area)


#change values of 999.99 to missing value (Perc_Rip and Perc_Pre1914)

ds[ds==999.99] <- NA

#extract records that are for farmland

ds = ds[ds$farm_flag==1,]

#check data
head(ds)

#check for missing values

any(is.na(ds$delta_tvp)) == F # delta_tvp does have NAs, but this function is not showing that, 
				#if you pull out just delta_tvp (tvp=ds$delta_tvp) and run this logical it will report true, 
				# however, lme4 seems to handle missing values automatically

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#START MODELING

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

# extract the deviance of the model
devcomp = getME(null,"devcomp")
dev = as.numeric(devcomp$cmp[8])
dev   # to compare nested models take difference of dev and difference in 
	# estimated parameters and check Chi_squared distribution for significance

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#TESTING LEVEL-2 PREDICTORS of INTEREST through direct and indirect effects
#Although devstat is calculated Pr value should be sufficent for testing fixed effects


#add wr density as predictor of int (level 2 predictor)
M1.model = lmer(delta_tvp ~ WR_density + (1 | HUC12), data = ds, REML = "False") 
summary(M1.model)
devcomp = getME(M1.model,"devcomp")
dev.M1 = as.numeric(devcomp$cmp[8])
dev.M1
devstat = dev -dev.M1
devstat #difference in number of est parameters (q) should be 1, this is not significant at 0.05 lvl

#add Rip as predictor of int (level 2 predictor)
M2.model = lmer(delta_tvp ~ (Rip) + (1 | HUC12), data = ds, REML = "False") 
summary(M2.model)
devcomp = getME(M2.model,"devcomp")
dev.M2 = as.numeric(devcomp$cmp[8])
dev.M2
devstat = dev -dev.M2
devstat #difference in number of est parameters (q) should be 2, model fit not improved but fixed effects sig

#add P1914  as predictor of int (level 2 predictor)
M3.model = lmer(delta_tvp ~ (P1914) + (1 | HUC12), data = ds, REML = "False") 
summary(M3.model)
devcomp = getME(M3.model,"devcomp")
dev.M3 = as.numeric(devcomp$cmp[8])
dev.M3
devstat = dev -dev.M3
devstat #difference in number of est parameters (q) should be 2, model fit not improved but fixed effects sig

#add RIP  as interaction with wr density as predictor of int (level 2 predictor)
M4.model = lmer(delta_tvp ~  (WR_density*Rip) + (1 | HUC12), data = ds, REML = "False") 
summary(M4.model)
devcomp = getME(M4.model,"devcomp")
dev.M4 = as.numeric(devcomp$cmp[8])
dev.M4
devstat = dev.M1 -dev.M4
devstat #difference in number of est parameters (q) should be 1, model fit not improved but fixed effects sig

#add gw density  as predictor of int (level 2 predictor)
M5.model = lmer(delta_tvp ~ GW_dnsty15 + (1 | HUC12), data = ds, REML = "False") 
summary(M5.model)
devcomp = getME(M5.model,"devcomp")
dev.M5 = as.numeric(devcomp$cmp[8])
dev.M5
devstat = dev -dev.M5
devstat #difference in number of est parameters (q) should be 1, model fit sig improved, most fixed effects still sig

#add Avg_WSEL_5yrChange to model above
M6.model = lmer(delta_tvp ~  GW_dnsty15 + Avg_WSEL_5yrChange +(1 | HUC12), data = ds, REML = "False") 
summary(M6.model)
devcomp = getME(M6.model,"devcomp")
dev.M6 = as.numeric(devcomp$cmp[8])
dev.M6
devstat = dev.M5 -dev.M6
devstat #difference in number of est parameters (q) should be 1, not sig

#add Avg_WSEL_5yrChange to model above as interaction
M7.model = lmer(delta_tvp ~  GW_dnsty15*Avg_WSEL_5yrChange +(1 | HUC12), data = ds, REML = "False") 
summary(M7.model)
devcomp = getME(M7.model,"devcomp")
dev.M7 = as.numeric(devcomp$cmp[8])
dev.M7
devstat = dev.M6 -dev.M7
devstat #difference in number of est parameters (q) should be 2,not sig

#MAX WSEL not sig either

#add WR_density back to model with GW
M8.model = lmer(delta_tvp ~  GW_dnsty15 + WR_density +(1 | HUC12), data = ds, REML = "False") 
summary(M8.model)
devcomp = getME(M8.model,"devcomp")
dev.M8 = as.numeric(devcomp$cmp[8])
dev.M8
devstat = dev.M5 -dev.M8
devstat #difference in number of est parameters (q) should be 2,not sig

#add Rip back to model with GW
M8.model = lmer(delta_tvp ~  GW_dnsty15 + Rip +(1 | HUC12), data = ds, REML = "False") 
summary(M8.model)
devcomp = getME(M8.model,"devcomp")
dev.M8 = as.numeric(devcomp$cmp[8])
dev.M8
devstat = dev.M5 -dev.M8
devstat #difference in number of est parameters (q) should be 2,not sig

#add P1914 back to model with GW
M8.model = lmer(delta_tvp ~  GW_dnsty15 + P1914 +(1 | HUC12), data = ds, REML = "False") 
summary(M8.model)
devcomp = getME(M8.model,"devcomp")
dev.M8 = as.numeric(devcomp$cmp[8])
dev.M8
devstat = dev.M5 -dev.M8
devstat #difference in number of est parameters (q) should be 2,not sig

#add Rip*WR_density back to model with GW
M9.model = lmer(delta_tvp ~  GW_dnsty15 + WR_density*Rip +(1 | HUC12), data = ds, REML = "False") 
summary(M9.model)
devcomp = getME(M9.model,"devcomp")
dev.M9 = as.numeric(devcomp$cmp[8])
dev.M9
devstat = dev.M5 -dev.M9
devstat #difference in number of est parameters (q) should be 3, sig at 0.025


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#ADDING CONTROLS to PRIMARY LEVEL-2 PREDICTOR
#Use devstat for testing significance of random effects

#add lu change as level-2 control to gw_density, wr type model above
M10.model = lmer(delta_tvp ~  GW_dnsty15 + WR_density*Rip  + delta_lc +(1 | HUC12), data = ds, REML = "False") 
summary(M10.model)
devcomp = getME(M10.model,"devcomp")
dev.M10 = as.numeric(devcomp$cmp[8])
dev.M10
devstat = dev.M9 -dev.M10
devstat #difference in number of est parameters (q) should be 1, not sig

#add lu change as level-1 control to gw_density, wr type model above
M11.model = lmer(delta_tvp ~  GW_dnsty15 + WR_density*Rip  + delta_lc +(1 + delta_lc | HUC12), data = ds, REML = "False") 
summary(M11.model)
devcomp = getME(M11.model,"devcomp")
dev.M11 = as.numeric(devcomp$cmp[8])
dev.M11
devstat = dev.M10 -dev.M11
devstat #difference in number of est parameters (q) should be 1, random effect is significant

#add diversity as additional lvl-2 control 
M12.model = lmer(delta_tvp ~ GW_dnsty15 + WR_density*Rip  + delta_lc + diverse +(1 + delta_lc | HUC12), data = ds, REML = "False") 
summary(M12.model)
devcomp = getME(M12.model,"devcomp")
dev.M12 = as.numeric(devcomp$cmp[8])
dev.M12
devstat = dev.M11 -dev.M12
devstat #difference in number of est parameters (q) should be 1, not sig

#add diversity as interaction control, use na.action = na.exclude for easing use of residual outputs
M13.model = lmer(delta_tvp ~ GW_dnsty15 + WR_density*Rip + delta_lc*diverse +(1 + delta_lc | HUC12), data = ds, REML = "False", na.action=na.exclude) 
summary(M13.model)
devcomp = getME(M13.model,"devcomp")
dev.M13 = as.numeric(devcomp$cmp[8])
dev.M13
devstat = dev.M11 -dev.M13
devstat #difference in number of est parameters (q) should be 2, significant model fit improvement

#wr have become non-sig, remove interaction first
M14.model = lmer(delta_tvp ~ GW_dnsty15 + WR_density +Rip + delta_lc*diverse +(1 + delta_lc | HUC12), data = ds, REML = "False") 
summary(M14.model)
devcomp = getME(M13.model,"devcomp")
dev.M14 = as.numeric(devcomp$cmp[8])
dev.M14
devstat = dev.M14 - dev.M13
devstat #difference in number of est parameters (q) should be 2, significant model fit improvement

#model fit not improved, and sig values for wr_density and Rip go way up so stick with model 13


#####BEST MODEL - MODEL 13 #######
#interpretation: Increased gw density has positive effect on HUC average delta_tvp (largest effect size). 
#Increased water rights density has a small positive effect on tvp.
#The count of Pre1914 wr in a HUC has a tiny positve effect on tvp and moderates the effect of wr density marginally, where higher counts lead to a slightly reduced effect of wr density on tvp.
#should not really be interpreting our controls
#A change in landuse has a small positive effect on tvp. The slope of tvp with lu change varies significantly between HUCs.
#Diversity moderates the effect of landuse change, where increased diversity reduces the effect of landuse change on tvp.

#r = resid(M13.model)

#r2 = ranef(M13.model)



#####Extract and Plot and Analyze Residuals #########
r2 = HLMresid(M13.model,level = "HUC12",  type = "EB") # raw group level residuals using empirical bayesinstall.packages
r1 = HLMresid(M13.model,level = 1)#raw level 1 residuals
r2=data.frame(r2)
ID = rownames(r2) #convert the index value to a rowname to get HUC12 ID in column
r2=mutate(r2,ID = as.numeric(ID))
fv = fitted.values(M13.model)

dat = data.frame(ds$lat, ds$lon, ds$HUC12, r1, fv)#combine lat, lon, and residuals
dat = left_join(dat,r2,by =c("ds.HUC12" = "ID")) #join Lvl-2 residuals based on HUC12 ID
newdat=na.omit(dat) #omit points with missing lat/lon


coordinates(newdat) = c('ds.lon','ds.lat')#set coordinate system
bubble(newdat,zcol='r1', maxsize = 1)#bubble chart of lvl-1 residuals for SJ
bubble(newdat,zcol='delta_lc1', maxsize = 1)#bubble chart of lvl-2 residuals (for delta_lc1 rand effects) for SJ 
bubble(newdat,zcol='X.Intercept.', maxsize = 1)#bubble chart of lvl-2 residuals (for intercept rand effects) for SJ
bubble(newdat,zcol='fv')#bubble chart of lvl-1 fitted values for SJ

var1.mod = variogram(r1~1,data=newdat,alpha=c(0,45,90,135), cutoff = 0.1) #variogram for level-1 residuals
plot(var1.mod)

var2.mod = variogram(X.Intercept.~1,data=newdat,alpha=c(0,45,90,135), cutoff = 0.5) #variogram for level-2 residuals
plot(var2.mod)

var2b.mod = variogram(delta_lc1~1,data=newdat,alpha=c(0,45,90,135), cutoff = 0.5) #variogram for level-2 residuals
plot(var2b.mod)

var2c.mod = variogram(delta_lc1~1,data=r2, cutoff = 0.5) #variogram for level-2 residuals, need to first define a coordinate system for the HUCs
plot(var2c.mod)

#moran's I and local moran's I
#use nlme cor structure to account for autocorr in lvl 1
