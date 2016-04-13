
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
setwd("C:/Users/tuan/Dropbox/Kate_Emily/Data")
#set up data, farming pixels only
dsf <- tbl_df(read.csv("pt_data_kn.txt", stringsAsFactors = FALSE))
ds = dsf[dsf$farm_flag==1,]
ds = transform(ds, 
               FID = as.numeric(FID),
               tvp_09 = as.numeric(tvp_09),
               tvp_14 = as.numeric(tvp_14),
               delta_tvp = as.numeric(delta_tvp),
               delta_lc = as.numeric(delta_lc),
               HUC12 = as.numeric(HUC12),
               AreaAcres = as.numeric(AreaAcres),
               GW_dnsty15 = as.numeric(GW_dnsty15),
               WR_density = as.numeric(WR_density))
ds[ds==999.99] <- NA
data <- ds[ds$farm_flag == 1,]
Rip <- ds$COUNT_Rip/ds$fmmp_area
P1914 <- ds$COUNT_Pre1914/ds$fmmp_area
y <- ds$tvp_14

#build group index that corresponds to adjacency matrix
uniq<-unique(ds$HUC12)
n.HUC <- length(uniq)
HUC <- rep(NA,dim(ds)[1])
for (i in 1:n.HUC){
  HUC[ds$HUC12 ==uniq[i]]<-i
}
HUC_index <- cbind(uniq, 1:56)

dsnew <- mutate(ds,HUCIDa=HUC, HUCID=HUC)

# adjacency matrix at group level
H <- inla.read.graph(filename = "shp.graph") 

#building a mesh for the level 1 data
coords<- select(dsnew,lon, lat) #pull coordinates
coordinates(coords) = c('lon','lat')
proj4string(coords)<-   #specify the projection
  " +init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83
+no_defs +towgs84=0,0,0"
coords<-as.data.frame(coords) #convert to workable formate
mesh0 <-inla.mesh.2d(coords, max.edge=0.1, offset=c(0.1,0.1), cutoff=0.01) 
plot(mesh0)
#for mesh1 added inner&outer max.edge with larger outer edge, outer triangles got bigger
mesh1 <-inla.mesh.2d(coords, max.edge=c(0.1,0.2), offset=c(0.1,0.1), cutoff=0.01) 
plot(mesh1)
#for mesh2 increase the cutoff to increase size of inner triangles
mesh2 <-inla.mesh.2d(coords, max.edge=c(0.1,0.2), offset=c(0.1,0.1), cutoff=0.02) 
plot(mesh2)
#for mesh3 remove the outer domain by removing second value in max.edge and offset lists
mesh3 <-inla.mesh.2d(coords, max.edge=c(0.1), offset=c(0.1), cutoff=0.02) 
plot(mesh3)

#creating a model fitting object from mesh (Matern SPDE object simulates values at the vertices)
spde <-inla.spde2.matern(mesh=mesh3, alpha=2) #pg 206, this mesh is for accounting for continuous spatial processes (landscape,etc..), and not for interpolating predictor values
# creating the observation/projection matrix (sparse weight matrix pg 205)
coords <-cbind(coords$lon, coords$lat)
A.est <-inla.spde.make.A(mesh=mesh3,loc=coords)

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

#multilevel model with land use factors e
formulae.mlm <-y ~ 1 + factor(lc_14) + f(HUC12, model="iid")
outpute.mlm <- inla (formulae.mlm, family = "gaussian", data=ds)
mlme.summary <- summary (outpute.mlm)

#multilevel model with land use factors f
formulaf.mlm <-y ~ 1 + GW_dnsty15*factor(lc_14) + f(HUC12, model="iid")
outputf.mlm <- inla (formulaf.mlm, family = "gaussian", data=ds)
mlmf.summary <- summary (outputf.mlm)

#multilevel model y
formulay.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc + f(HUC12, model="iid")+ f(HUC12a, FID, model="iid", constr=TRUE)
outputy.mlm <- inla (formulay.mlm, family = "gaussian", data=dsnew)
mlmy.summary <- summary (outputy.mlm)

#spatial model no predictors
formula.spat <-y ~ f(HUCID, model="bym",graph=H, scale.model=TRUE) #page 183 for scale and pg 182 for auto priors with bym
output.spat <- inla (formula.spat, family ="gaussian", data=dsnew)
mod.spat <-summary(output.spat)
head(output.spat$summary.random$HUCID)  #this gives us random effects for HUCs (0:n.HUC area-specific residuals, n.HUC:2*n.HUC are spatially structured residuals)
head(output.spat$summary.fixed) #fixed effects, same as summary, grand mean

#multilevel model z, this where we add an iid effect to the CAR besag model (same as bym model with no predictors)
formulaz.mlm <-y ~  f(HUCID, model="besag", graph = H) + f(HUCIDa, model="iid", constr=TRUE)
outputz.mlm <- inla (formulaz.mlm,  data=dsnew)
mlmz.summary <- summary (outputz.mlm) #high precision mean on hyperparameters may be due to few neighboring hucs (Pg 178)

#multilevel model x (with predictors)
formulax.mlm <-y ~ 1 + GW_dnsty15 + WR_density + f(HUCID, model="besag", graph = H)
outputx.mlm <- inla (formulax.mlm,  data=dsnew)
mlmx.summary <- summary (outputx.mlm)

#multilevel model j (with predictors), this where we add an iid effect to the CAR besag model, same as bym model with predictors
formulaj.mlm <-y ~ 1 + GW_dnsty15 + WR_density + f(HUCID, model="besag", graph = H) + f(HUCIDa, model="iid", constr=TRUE)
outputj.mlm <- inla (formulaj.mlm,  data=dsnew)
mlmj.summary <- summary (outputj.mlm)

#multilevel model k (with predictors), random and spatially structured effects
formulak.mlm <-y ~ 1 + GW_dnsty15 + WR_density + f(HUCID, model="bym", graph = H, scale.model=TRUE)
outputk2.mlm <- inla (formulak.mlm,  data=dsnew)
mlmk2.summary <- summary (outputk2.mlm)

#multilevel model m (with predictors), random and spatially structured effects, with lc factors
formulam.mlm <-y ~ 1 + GW_dnsty15 + WR_density + factor(lc_14) + f(HUCID, model="bym", graph = H, scale.model=TRUE)
outputm.mlm <- inla (formulam.mlm,  data=dsnew)
mlmm.summary <- summary (outputm.mlm)

#multilevel model t, this where we add second indexing , IGNORE THIS ONE I THINK...
formulat.mlm <-y ~ 1 + GW_dnsty15 + WR_density + f(HUCID, model="besag", graph = H) + f(HUCIDa, FID, model="iid", constr=TRUE)
outputt.mlm <- inla (formulat.mlm,  data=dsnew)
mlmt.summary <- summary (outputt.mlm)


##page 206, model fitting with mesh
##use tvp_14 as y
##add lc_14 categories
##add random slope, pg 160-161, f(id2,x,model="iid")

#3-9-2016 just a spatial model using mesh
formula.mlm <-y ~ -1 + intercept + f(spatial.field, model = spde)
data = list(y=y, intercept=rep(1,spde$n.spde), spatial.field=1:spde$n.spde) 
output.mlm <- inla (formula.mlm,  data=data, control.predictor =list(A=A.est, compute=TRUE))
mlm.summary <- summary (output.mlm)


#3-9-2016 updating multilevel model m with mesh, for simplicity use delta_lc instead of factors
formula.mlm <-y ~ 1 + GW_dnsty15 + WR_density + delta_lc + f(HUCID, model="bym", graph = H, scale.model=TRUE) + f(spatial.field, model = spde)
data = list(y=y, GW_dnsty15= dsnew$GW_dnsty15, WR_density=dsnew$WR_density, delta_lc=dsnew$delta_lc, HUCID=dsnew$HUCID, spatial.field=1:spde$n.spde) 
output.mlm <- inla (formula.mlm,  data=data, control.predictor =list(A=A.est, compute=TRUE))
mlm.summary <- summary (output.mlm)