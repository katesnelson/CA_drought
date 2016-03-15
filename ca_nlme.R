
library(memisc)
library(dplyr)
library (lme4)
library (lmerTest)
library(HLMdiag)
library (gstat)
library(sp)
library(nlme)
library(spdep)
library(rgdal)

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
devstat 
#model fit not improved, and sig values for wr_density and Rip go way up so stick with model 13


#####BEST MODEL - MODEL 13 #######
#interpretation: Increased gw density has positive effect on HUC average delta_tvp (largest effect size). 
#Increased water rights density has a small positive effect on tvp.
#The count of Pre1914 wr in a HUC has a tiny positve effect on tvp and moderates the effect of wr density marginally, where higher counts lead to a slightly reduced effect of wr density on tvp.
#should not really be interpreting our controls
#A change in landuse has a small positive effect on tvp. The slope of tvp with lu change varies significantly between HUCs.
#Diversity moderates the effect of landuse change, where increased diversity reduces the effect of landuse change on tvp.



#####Extract and Plot and Analyze Residuals #########
r2 = HLMresid(M13.model,level = "HUC12",  type = "EB") # raw group level residuals using empirical bayesinstall.packages
r = HLMresid(M13.model,level = 1)#raw level 1 residuals
r2=data.frame(r2)
ID = rownames(r2) #convert the index value to a rowname to get HUC12 ID in column
r2=mutate(r2,ID = as.numeric(ID))
fv = fitted.values(M13.model)

dat = data.frame(ds$lat, ds$lon, ds$HUC12, r)#combine lat, lon, and residuals
dat = left_join(dat,r2,by =c("ds.HUC12" = "ID")) #join Lvl-2 residuals based on HUC12 ID
newdat=na.omit(dat) #omit points with missing lat/lon


coordinates(newdat) = c('ds.lon','ds.lat')#set coordinate system
bubble(newdat,zcol='r', maxsize = 1)#bubble chart of lvl-1 residuals for SJ
bubble(newdat,zcol='delta_lc1', maxsize = 1)#bubble chart of lvl-2 residuals (for delta_lc1 rand effects) for SJ 
bubble(newdat,zcol='X.Intercept.', maxsize = 1)#bubble chart of lvl-2 residuals (for intercept rand effects) for SJ
bubble(newdat,zcol='fv')#bubble chart of lvl-1 fitted values for SJ
proj4string(newdat)<- 
" +init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83
+no_defs +towgs84=0,0,0"

var1.mod = variogram(r~1,data=newdat,alpha=c(0,45,90,135)) #variogram for level-1 residuals
plot(var1.mod)

var2.mod = variogram(r~1,data=newdat,alpha=c(0,45,90,135), width=.55) #variogram for level-2 residuals
plot(var2.mod)

var2b.mod = variogram(delta_lc1~1,data=newdat,alpha=c(0,45,90,135)) #variogram for level-2 residuals
plot(var2b.mod)

var2c.mod = variogram(X.Intercept.~1,data=newdat,alpha=c(0,45,90,135)) #variogram for level-2 residuals, need to first define a coordinate system for the HUCs
plot(var2c.mod)

#moran's I and local moran's I
#use nlme cor structure to account for autocorr in lvl 1

#############Using NLME################
ca.model <- lme (fixed = delta_tvp ~ GW_dnsty15 + WR_density*Rip + delta_lc*diverse, random =  ~ 1 + delta_lc | HUC12, method = "ML", data = ds, na.action = na.omit)
summary (ca.model) #model estimates same as what we got using lmer

#wow! the model is 1.1Gb in size, almost no difference between estimates from M13 model, and REALLY SMALL phi
ca.gau <- update (ca.model, correlation = corGaus(1, form = ~ lat + lon | HUC12), method = "ML",control=lmeControl(opt = "optim") )
summary(ca.gau)

#also 1.1Gb, estimates are signifcantly different!!, does a good job of smoothing out residual values within HUCs, still sig differences between HUCS
ca.exp <- update (ca.model, correlation = corExp(1, form = ~ lat + lon | HUC12), method = "ML",control=lmeControl(opt = "optim") )
summary(ca.exp)

# same as ca.exp, automatically groups the residual autocorr in HUC12 (probably because the random effects are)
ca.expb <- update (ca.model, correlation = corExp(1, form = ~ lat + lon), method = "ML",control=lmeControl(opt = "optim") )
summary(ca.expb)

#####Extract and Plot NLME Residuals ###########
r1 = resid(ca.expb)
r2 = ranef(ca.exp)
fv = fitted.values(ca.exp)

head(r1)
head(r2)

dsnew = select_(ds, 'lat', 'lon', 'HUC12', 'delta_tvp', 'GW_dnsty15','WR_density','Rip','delta_lc','diverse')
dsnewb = na.omit(dsnew)

r2=data.frame(r2)
ID = rownames(r2) #convert the index value to a rowname to get HUC12 ID in column
r2=mutate(r2,HUC12 = as.numeric(ID))


dat.exp = data.frame(dsnewb$lat, dsnewb$lon, dsnewb$HUC12, r1)#combine lat, lon, and residuals
da.expt = left_join(dat,r2,by =c("dsnewb.HUC12" = "HUC12")) #join Lvl-2 residuals based on HUC12 ID
newdat.exp=na.omit(dat.exp) #omit points with missing lat/lon


coordinates(newdat.exp) = c('dsnewb.lon','dsnewb.lat') #set coordinate system
proj4string(newdat.exp)<- 
  " +init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83
+no_defs +towgs84=0,0,0"
bubble(newdat.exp,zcol='r1', maxsize = 1)#bubble chart of lvl-1 residuals for SJ
varexpb.mod = variogram(r1~1,data=newdat.exp,alpha=c(0,45,90,135), width=.01) #variogram for level-2 residuals
plot(varexpb.mod)
