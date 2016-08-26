

##################BEST MODELS SO FAR############################

sampleVar <- variogramST(ELEV~1, sto, tunit="days", tlags=seq(from=0,to=180, by=15), cutoff = 60, na.omit=T)  #shorter time range to address non-stationarity through time, space already truncated to address same
#write.table(sampleVar, 'sampleVar.csv', sep=",")
#sampleVar <- read.table('sampleVar.csv', sep=",")

#diagnostic plots
acf(sto@data, 'xts')  #up to 200 days significant
plot(sampleVar, map=F) #range of around 60, nugget less than 1
plot(sampleVar$timelag, sampleVar$gamma, map=F )
plot(sampleVar) #range of around 60, nugget less than 1
plot(sampleVar, wireframe=T) #visualize space/time lag
#for unprojected great circle distances (km).
#gamma values conditional on treatment of data (log, etc)

#estimate space-time anisotropy, or directional dependence of data
#http://finzi.psych.upenn.edu/library/gstat/html/estiStAni.html
t<-estiStAni(sampleVar, interval = c(1,85), method = "metric", vgm(1000000, "Gau", 60, 13000), s.range = 60, t.range = 60)
#data too sparse?  http://stats.stackexchange.com/questions/35316/problems-estimating-anisotropy-parameters-for-a-spatial-model

#parameter range
#paramcheck:http://giv-graeler.uni-muenster.de:3838/spacetime/
pars.l <- c(sill.s = 0.1, range.s = 10, nugget.s = 0.1, sill.t = 0.1, range.t = 1, 
            nugget.t = 0.1, sill.st = 0.1, range.st = 10, nugget.st = 0.1, anis = 0.1)
pars.u <- c(sill.s = 2000000, range.s = 60, nugget.s = 200000, sill.t = 500000, range.t = 180, 
            nugget.t = 200000, sill.st = 3000000, range.st = 100000, nugget.st = 100000, anis = 1000) #upper anis based on estiStAni

##################
#metric model
##################
#http://geostat-course.org/system/files/part01.pdf
#Temporal domain rescaled to match spatial one.  Lag classes include space, time, and space-time distances.  
#Using a single model forces all distances to behave in the same way (unlikely)
#assumes identical covariance unctions for space and time
#includes space-time anisitropy that allows some flexibility
#all distances treated equally, only one joint variogram for all three.

metVar <- vgmST("metric",
                joint = vgm(1000000,"Mat",80,114000, kappa=10), #sill, model, range, nugget
                stAni = 180) #directional dependence
metVarDF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=0)
MSEmetVarDF <- attr(metVarDF, "MSE")
metVarBF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u, stAni = 180)
MSEmetVarBF <- attr(metVarBF, "MSE") 

#fit.model = 7 is method that controls the weighing of the residuals between empirical and model surface;
#7 is the ST analog to commonly used spatial weighting; higher confidence in distances with many pairs

###################
#product-sum
##################
#relax assummptions in covariance matrix, does not assume separability
#need to provude both spatial and temporal components, plus k, weighting product

psVar <- vgmST("productSum", 
               space=vgm(180, "Gau", 60, 3),   #sill, model, range, nugget
               time=vgm(70, "Sph", 60, 1),
               k=700)
psVarDF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=0)  #data fit
MSEpsVarDF <- attr(psVarDF, "MSE")
psVarBF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=7, control = list(parscale = c(1,10,1,1,0.1,1,10)),
                           lower = rep(0.0001, 70), stAni = 180)
#k must be positive error
MSEpsVarBF <- attr(psVarBF, "MSE")

####################
#sum-metric
####################
#metric with spatial and temporal covariance models plus joint component with anisitropy
#maximum flexibility

smVar <- vgmST("sumMetric", 
               space=vgm(1000000, "Gau", 50, 20000),   #sill, model, range, nugget
               time=vgm(91000, "Sph", 50, 13000),
               joint=vgm(13000,"Sph",20000, 22000),
               stAni=180)
smVarDF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=0)  
MSEsmVarDF <- attr(smVarDF, "MSE")
smVarBF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u,control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                                                                                                        maxit=1e4), stAni =180)
MSEsmVarBF <- attr(smVarBF, "MSE")

######################
#simple sum-metric
######################
#like sum metric but restrict all components to have a single nugget

ssmVar <- vgmST("simpleSumMetric", 
                space=vgm(1000000, "Gau", 60),   #sill, model, range, nugget
                time=vgm(70000, "Sph", 42),
                joint=vgm(1,"Sph",80),
                nugget=74000, stAni=180)
ssmVarDF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=0)  
MSEssmVarDF <- attr(ssmVarDF, "MSE")
ssmVarBF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u,control = list(parscale = c(1,10,1,1,1,100,1,10)), stAni=180)
MSEssmVarBF <- attr(ssmVarBF, "MSE")

#####################
#separable model
#####################

#not using since bad assumptions, but documentation and scripts below:
#http://geostat-course.org/system/files/part01.pdf: assumptions of isotropy and stationarity
#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/: pararmeter selection for vgmST doesn't seem to matter so much 
# sepVar <- vgmST("separable", 
#                   space=vgm(15, "Exp", 60, 0),   #sill, model, range, nugget
#                   time=vgm(15, "Exp", 40, 0),
#                   sill = 30,
#                   stAni = 100)
# 
# #check how well the model fits data with fit.method=0 which keeps the model but calculates MSE compared to actual data
# sepVarDF <- fit.StVariogram(sampleVar, sepVar, method="L-BFGS-B", fit.method=0)  #data fit
# MSEsepVarDF <- attr(sepVarDF, "MSE")  #think error of the eye fit
# sepVarBF <- fit.StVariogram(sampleVar, sepVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
# MSEsepVarBF <- attr(sepVarBF, "MSE") 
#extractPar(sepVarBF)


##########################################
#MODEL COMPARISON
###########################################

#space-time-kriging document
#sampleVar
#plot empirical vars
#create version of table 2, figure 4 shows space-time biases
#find model with smallest RMSE

#plot all 2D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF, psVarBF), all=T) 

#plot all 3D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF, psVarBF), all=T, wireframe=T, zlim=c(0,1000000),  
     zlab=NULL, xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)))

#plot all diff 3D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF, psVarBF), all=T, wireframe=T, zlim=c(-100000,1000000),  
     zlab=NULL, xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)), diff=TRUE)

#plot fit metrics (data fit)
MSE_DF <- c(MSEmetVarDF, MSEsmVarDF, MSEssmVarDF, MSEpsVarDF) 
barplot(MSE_DF, main="MSE of data fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Sum-metric", "Product-sum"))

#plot fit metrics (best fit)
MSE_BF <- c(MSEmetVarBF, MSEsmVarBF, MSEssmVarBF, MSEpsVarBF)  
barplot(MSE_BF, main="MSE of best fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Sum-metric", "Product-sum"))

extractPar(ssmVarBF) #gives fitted parameter values, anis is always the max?
extractPar(smVarBF)
extractPar(metVarBF)
extractPar(psVarBF)


###For kriging and plotting#######




#create new spatiotemporal prediction grid

#create spatial reference
newproj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
rast <- raster()
extent(rast) <- extent(d)
rast <- as(rast, 'SpatialPixels') #360 x 180 cells
res(rast@grid@cells.dim) <- c(180,90)

###add line here to clip grid to boundaries

#create time object
#tm.grid <- dfs$TIME[order(dfs$TIME)]
tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 12*15)

tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 3)
#sto frame
grid.st <- STF(rast, tm.grid)

sto@sp=spTransform(sto@sp,'+init=epsg:3310')  # transform projection for kriging to NAD83/california albers
grid.st=spTransform(grid.st,sto@sp@proj4string) #set to same proj as sto

##Interpolation
#nmax of 20 yields very large (spatial scale) trends that do not match the resolution we want well, nmax of 10 seems better
#stAni: taking anis value directly from model parameters gives really bad results, Edzer P does this "stAni = fitMetricModel$stAni/24/3600" 
            #where the units of stani are km/day, yielding an stani in km/sec (not sure why they do this so assume it may be a standard function setting)
pred<-krigeST(ELEV~1, data=sto, newdata=grid.st, nmax = 10, modelList=smVarBF, stAni = smVarBF$stAni/24/3600, progress=TRUE) 

## Plot predictions
coordinates(dfs) <- ~LON+LAT
dfs@proj4string <- CRS('+init=epsg:4326')
dfs=spTransform(dfs,sto@sp@proj4string)  #prepare observed values for overlay on prediction grid
layout <-list("sp.points", dfs,  first=F, cex=.1)
stplot(pred, sp.layout = layout)

#Pull pred values for comparison
b<-buffer(dfs,5000) #buffer the pts (5km) for overlay as pred3 spatial object is listed as pts too
#if we can convert the pred to a spatial grid this might work better, not sure what an optimal buffer distance would be

select <- pred3[b,1,"var1.pred"] #select prediction values that overlay with points (for first time instance)
#will want to write something such that we look at a set of obs for a specific time frame and use that same 
#frame for the prediction pull