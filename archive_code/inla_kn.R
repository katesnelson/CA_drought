
library(memisc)
library(dplyr)
library (lme4)
library (lmerTest)
library(HLMdiag)
library (gstat)
library(sp)
library(nlme)
library(rgdal)
library(INLA)
library(spdep)
setwd("C:/Users/nelsonks/Dropbox/Kate_Emily/Data")
#set up data, farming pixels only
dsf <- tbl_df(read.csv("pt_data_kn.txt", stringsAsFactors = FALSE))
ds = dsf[dsf$farm_flag==1,]
ds = transform(ds, 
               FID = as.character(FID),
               tvp_09 = as.numeric(tvp_09),
               tvp_14 = as.numeric(tvp_14),
               delta_tvp = as.numeric(delta_tvp),
               delta_lc = as.numeric(delta_lc),
               HUC12 = as.character(HUC12),
               AreaAcres = as.numeric(AreaAcres))
ds[ds==999.99] <- NA
data <- ds[ds$farm_flag == 1,]
Rip <- ds$COUNT_Rip/ds$fmmp_area
P1914 <- ds$COUNT_Pre1914/ds$fmmp_area
y <- ds$delta_tvp


#pooled/aggregated data model
pooled <- y ~ 1 + GW_dnsty15 + WR_density + delta_lc
output.pooled <- inla(pooled, family = "gaussian", data = ds, verbose=TRUE) #status 253 error
pool.summary <- output.pooled$summary.fixed

#independent unit model (ind. model for each HUC)
ind <- y ~ 1 + factor(HUC12)+ GW_dnsty15 + WR_density + delta_lc 
output.ind <- inla(ind, family = "gaussian", data = ds, verbose=TRUE) 
ind.summary <- output.ind$summary.fixed

#multilevel model
formula.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc + f(HUC12, model="iid")
output.mlm <- inla (formula.mlm, family = "gaussian", data=ds)
mlm.summary <- summary (output.mlm)

#multilevel model a
formulaa.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc*diverse + f(HUC12, model="iid")
outputa.mlm <- inla (formulaa.mlm, family = "gaussian", data=ds)
mlma.summary <- summary (outputa.mlm)

#multilevel model b
formulab.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc + f(HUC12, model="iid", hyper=list(prec=list(prior="loggamma",param=c(1,0.001))))
outputb.mlm <- inla (formulab.mlm, family = "gaussian", data=ds, control.predictor=list(compute=TRUE),control.fixed=list(mean=0, prec=0.0001, mean.intercept=0, prec.intercept=0.0001))
mlmb.summary <- summary (outputb.mlm)

#multilevel model c
formulac.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc + f(HUC12, model="iid", hyper=list(prec=list(prior="loggamma",param=c(0.1,0.001))))
outputc.mlm <- inla (formulac.mlm, family = "gaussian", data=ds, control.predictor=list(compute=TRUE),control.fixed=list(mean=0, prec=0.0001, mean.intercept=0, prec.intercept=0.0001))
mlmc.summary <- summary (outputc.mlm)

#multilevel model d
formulad.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc + f(HUC12, model="iid", hyper=list(prec=list(prior="loggamma",param=c(0.001,0.0001))))
outputd.mlm <- inla (formulad.mlm, family = "gaussian", data=ds, control.predictor=list(compute=TRUE),control.fixed=list(mean=0, prec=0.0001, mean.intercept=0, prec.intercept=0.0001))
mlmd.summary <- summary (outputd.mlm)

#multilevel model with land use factors
formulae.mlm <-y ~ 1 + GW_dnsty15 + WR_density + factor(lc_14) + f(HUC12, model="iid")
outpute.mlm <- inla (formulae.mlm, family = "gaussian", data=ds)
mlme.summary <- summary (outpute.mlm)

#multilevel model with land use factors
formulaf.mlm <-y ~ 1 + GW_dnsty15*factor(lc_14) + f(HUC12, model="iid")
outputf.mlm <- inla (formulaf.mlm, family = "gaussian", data=ds)
mlmf.summary <- summary (outputf.mlm)