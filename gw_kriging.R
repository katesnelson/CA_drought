library(dplyr)
library(gstat)
library(spacetime) 
library(sp)  
library(xts) 
library(rgdal)
library(raster)
library(Metrics)
library(ggplot2)
library(lattice)


#This script reads in groundwater elevation observations for the state of California extracted and cleaned in the script GroundwaterDataClean.R, builds a
#variogram of a subset of observed data (training data), fits a dependency model to the variogram, krigs data to a 10km x 10km, monthly interval grid using 
#the variogram models and training data, and produces diagnostic statistics and plots of the predicted groundwater elevations against held-out observed data.
#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/
#http://localhost/rstudio//files/R/x86_64-pc-linux-gnu-library/3.2/gstat/doc/st.pdf

#set.seed(10010)

#load data 
d_full <- as.data.frame(read.csv(paste('/home/kate/CA/data/gw.csv', sep=','), header=T), stringsAsFactors=F)
coordinates(d_full) <- ~LON+LAT
d_full@proj4string <- CRS('+init=epsg:4326') #confirm original projection

#subset point dataset with CV shapefile
d_full <- spTransform(d_full, CRS('+init=epsg:3310')) #change projection for clipping
cv <- shapefile('/home/kate/CA/data/cv_huc.shp')
cv <- spTransform(cv, CRS(proj4string(d_full))) #change projection for clipping
d <- d_full[cv,]

#create final dataframe
df <- data.frame(d@data$ELEV, as.Date(d@data$DATE), as.data.frame(d@coords[,1]), as.data.frame(d@coords[,2]))
names(df) <- c("ELEV", "TIME", "LON", "LAT")

#clean up outliers
df <- subset(df, df$ELEV >= -50 & df$ELEV < 4000, select=c("ELEV", "TIME", "LON", "LAT"))
constant_transform <- abs(min(df$ELEV)) + 1 #linear transformation to remove negative values
df$ELEV <- log(df$ELEV + constant_transform) #normally distributed

#select random subset of data (60%), does not necessarily cover all space-time points with observations
dfs <- df[sample(nrow(df), 0.6*nrow(df)), ] 

#remove any duplicates
dfs_clean <- dfs[!(duplicated(dfs[c("LAT", "LON", "TIME")]) | duplicated(dfs[c("LAT","LON", "TIME")], fromLast = TRUE)), ]

#sample of held out datapoints to test fit later
held_out <- rbind(df, dfs) #duplicate sampled observations
held_out <- held_out[! duplicated(held_out, fromLast=TRUE) & seq(nrow(held_out)) <= nrow(df), ]

#space-time object
sto <- stConstruct(dfs_clean, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval
sto@sp@proj4string <- CRS('+init=epsg:3310')  # set projection to original NAD83/california albers
sto <- spTransform(sto, CRS('+init=epsg:4326'))  #for the construction of the sampleVar, need to have unprojected (WGS84)
#stplot(sto)

#create sample variogram
#sampleVar <- variogramST(ELEV~1, sto, tunit="days", tlags=seq(from=0,to=180,by=15), cutoff = 60, na.omit=T, progress=T)  
#saveRDS(sampleVar, "sampleVar_0816.rds")
#test <- readRDS("sampleVar_0816.rds")

sampleVar2 <- readRDS("/home/kate/CA/sampleVar_0816.rds")

#clean up unnecessary files
remove(d_full, df, dfs)


###########################################################################################################################
#diagnostic plots
###########################################################################################################################
sampleVar<-sampleVar2
plot(sampleVar, map=F) 
plot(sampleVar, wireframe=T) 

#at spatial lag near 0, temporal autocorrelation (near zero, after 10 days)
timeVar <- sampleVar[sampleVar$spacelag == 2,]
plot(timeVar$timelag, timeVar$gamma)
acf(sto@data$ELEV, type="partial", lag.max = 100) 

#at time lag of zero, spatial autocorrelation
spVar <- sampleVar[sampleVar$timelag == 7.5,]
plot(spVar, map=F)

#contour plot
contourplot(gamma~spacelag+timelag, sampleVar, cuts=25)
#t<-estiStAni(sampleVar, interval = c(1,100), method = "metric", vgm(1000000, "Gau", 60, 13000), s.range = 60, t.range = 150)
#print(t)


############################################################################################################
#PARAMETER RANGE FOR ESTIMATION
############################################################################################################

#parameter range, http://giv-graeler.uni-muenster.de:3838/spacetime/, provide huge range here
pars.l <- c(sill.s = 0.1, range.s = 0.01, nugget.s = 0.01, sill.t = 0.1, range.t = 0.01, 
            nugget.t = 0.05, sill.st = 0.1, range.st = 1, nugget.st = 0.05, anis = 0.01)
pars.u <- c(sill.s = 1, range.s = 60, nugget.s = 0.4, sill.t = 1, range.t = 200, 
            nugget.t = 0.4, sill.st = 1, range.st = 180, nugget.st = 1, anis = 10)

############################################################################################################
#METRIC MODEL
############################################################################################################

#http://geostat-course.org/system/files/part01.pdf
#Temporal domain rescaled to match spatial one.  
#Lag classes include space, time, and space-time distances.  
#Using a single model forces all distances to behave in the same way (unlikely)
#assumes identical covariance functions for space and time
#includes space-time anisitropy that allows some flexibility
#all distances treated equally, only one joint variogram for all three.

#fit.model = 7 is method that controls the weighing of the residuals between empirical and model surface;
#7 is the ST analog to commonly used spatial weighting; higher confidence in distances with many pairs

# metVar <- vgmST("metric",
#                 joint = vgm(1,"Gau",80,0.2, kappa=0.6), #sill, model, range, nugget
#                 stAni = 0.01) #directional dependence
# #metVarDF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=0) #data fit
# #MSEmetVarDF <- attr(metVarDF, "optim.output")$value  #not working
# metVarBF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u,
#                             control = list(parscale = c(10, 20, 5, 10)))
# MSEmetVarBF <- attr(metVarBF, "optim.output")$value
# print(MSEmetVarBF)

############################################################################################################
#PRODUCT SUM
############################################################################################################

#relax assummptions in covariance matrix, does not assume separability
#need to provide both spatial and temporal components, plus k, weighting product

psVar <- vgmST("productSum", 
               space=vgm(.1, "Gau", 60, 0.2), 
               time=vgm(.1, "Exp", 80, 0.2),
               k=100)
#psVarDF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=0)  
#MSEpsVarDF <- attr(psVarDF, "optim.output")$value
psVarBF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u,
                           control = list(parscale = c(1, 10, 1, 1, .1, 1, 10), maxit=1e4), stAni = .01)
MSEpsVarBF <- attr(psVarBF, "optim.output")$value
print(MSEpsVarBF)

############################################################################################################
#SUM METRIC
############################################################################################################

#metric with spatial and temporal covariance models plus joint component with anisitropy
#maximum flexibility

smVar <- vgmST("sumMetric", 
               space=vgm(.1, "Gau", 60, 0.2), #sill, model, range, nugget
               time=vgm(.1, "Gau", 80, 0.2),
               joint=vgm(.1,"Exp",100,0.2),
               stAni=.01)
#smVarDF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=0)  
#MSEsmVarDF <- attr(smVarDF, "optim.output")$value
smVarBF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u, 
                           control = list(parscale = c(1, 100, 1, 1, 0.5, 1, 1, 100, 1, 100), maxit=1e4), stAni = .01)
MSEsmVarBF <- attr(smVarBF, "optim.output")$value
print(MSEsmVarBF)


############################################################################################################
#SIMPLE SUM METRIC
############################################################################################################

#like sum metric but restrict all components to have a single nugget

ssmVar <- vgmST("simpleSumMetric", 
                space=vgm(.1, "Gau", 60, 0.2),
                time=vgm(.1, "Gau", 80, 0.2),
                joint=vgm(.1,"Sph",100,0.2),
                nugget=.4, stAni=.01)
#ssmVarDF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=0)  
#MSEssmVarDF <- attr(ssmVarBF, "optim.output")$value
ssmVarBF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u, 
                            control = list(parscale = c(1, 10, 1, 1, 1, 100, 1, 10)), stAni=.01)
MSEssmVarBF <- attr(ssmVarBF, "optim.output")$value
print(MSEssmVarBF)

############################################################################################################
#MODEL COMPARISON
############################################################################################################

#plot all 2D
plot(sampleVar, list(smVarBF, ssmVarBF, psVarBF), all=T) 

#plot all 3D
plot(sampleVar, list(smVarBF, ssmVarBF, psVarBF), all=T, wireframe=T, 
     zlab=NULL, xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)))

#plot 3D differences
plot(sampleVar, list(smVarBF, ssmVarBF, psVarBF), all=T, wireframe=T, 
     zlab=NULL, xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)), diff=T)

#plot fit metrics (best fit)
#weighted MSE
MSE_BF <- c(MSEsmVarBF, MSEssmVarBF, MSEpsVarBF)  
barplot(MSE_BF, main="Weighted MSE", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Sum-metric", "Simple sum-metric", "Product-sum"))

#unweighted MSE
MSEssm <- attr(ssmVarBF, "MSE")
MSEsm <- attr(smVarBF, "MSE")
MSEps<- attr(psVarBF, "MSE")
MSE_BF_uw <- c(MSEsm, MSEssm, MSEps)
barplot(MSE_BF_uw, main="Unweighted MSE", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Sum-metric", "Simple sum-metric", "Product-sum"))


############################################################################################################
#PREDICTION GRID CONSTRUCTION
############################################################################################################

#WGS84 did not work, use projection

#create empty spatial grid
#coords in meters, 3310 projection
grd <- SpatialPixels(SpatialPoints(makegrid(cv, n=2900)), proj4string = proj4string(cv))
grd <- grd[cv,] #subset to central valley
#check resolution
#test_raster <- raster(grd) #10K x 10K
#print(test_raster)

#establish correct project for kriging
sto_krig <- spTransform(sto, CRS('+init=epsg:3310')) #project to krig

#create monthly time object
tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 12*15)

#empty sto object
grid.st <- STF(sp = as(grd, "SpatialPixels"), time = tm.grid)

#############################
#sum-metric
#############################
pred <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=20, modelList=smVarBF, stAni=smVarBF$stAni, progress=T) 
gridded(pred@sp) <- TRUE 

# NOTES on Spatial-Temporal Anisotropy from  Pebesma, 2015, gStat documentation. 
# In vgmSt stated that" A spatio-temporal anisotropy is the number of space units equivalent to one time unit." 
# Note that sumMetric, simpleSumMetric, and metric vgmST models require StAni as an input parameter, seperable and productSum models do not. 
# In krigeST stated that stAni is "a spatio-temporal anisotropy scaling assuming a metric spatio-temporal space. 
# Used only for the selection of the closest neighbours. This scaling needs only to be provided in case the model does not have a stAni parameter, 
# or if a different one should be used for the neighbourhood selection. Mind the correct spatial unit. Currently, no coordinate conversion is made 
# for the neighbourhood selection (i.e. Lat and Lon require a spatio-temporal anisotropy scaling in degrees per second)."

#adjustment to stAni 
pred2 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=20, modelList=smVarBF, stAni=smVarBF$stAni/24/3600, progress=T) #time unit of stAni appears to be annual by default (see gStat documentaiton) and varigram is in days
gridded(pred2@sp) <- TRUE 

#change spatial range, joint range and stAni 
smVarBF_t <- smVarBF
smVarBF_t$space$range <- smVarBF_t$space$range*1000 #conversion accounts for difference between units of variogram and model in WGS unprojected (km) and grid in projected (m)
smVarBF_t$joint$range <- smVarBF_t$joint$range*1000
smVarBF_t$stAni <- smVarBF_t$stAni*1000

pred3 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=20, modelList=smVarBF_t, stAni=smVarBF_t$stAni/24/3600, progress=T) 
gridded(pred3@sp) <- TRUE 

pred3.5 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=40, modelList=smVarBF_t, stAni=smVarBF_t$stAni/24/3600, progress=T) 
gridded(pred3@sp) <- TRUE 

#############################
#simple sum-metric
#############################

pred4 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=20, modelList=ssmVarBF, stAni=ssmVarBF_t$stAni, progress=T)
gridded(pred4@sp) <- TRUE

#adjustment to stAni
pred5 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=20, modelList=ssmVarBF, stAni=ssmVarBF$stAni/24/3600, progress=T)
gridded(pred5@sp) <- TRUE


#change spatial range, joint range and stAni
ssmVarBF_t <- ssmVarBF  #MAKE COPY
ssmVarBF_t$space$range <- ssmVarBF_t$space$range*1000
ssmVarBF_t$joint$range <- ssmVarBF_t$joint$range*1000
ssmVarBF_t$stAni <- ssmVarBF_t$stAni*1000

pred6 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=20, modelList=ssmVarBF_t, stAni=ssmVarBF_t$stAni/24/3600, progress=T) 
gridded(pred6@sp) <- TRUE 

pred6.5 <- krigeST(ELEV~1, data=sto_krig, newdata=grid.st, nmax=40, modelList=ssmVarBF_t, stAni=ssmVarBF_t$stAni/24/3600, progress=T) 
gridded(pred6@sp) <- TRUE 


###################################################
#compute res/rmse
###################################################
#build small residual datasets for model comparison

n <- 1000
pred_rmse <- krig_rmse(pred, n)
pred2_rmse <- krig_rmse(pred2, n)
pred3_rmse <- krig_rmse(pred3, n)
pred3.5_rmse <- krig_rmse(pred3.5, n)
pred4_rmse <- krig_rmse(pred4, n)
pred5_rmse <- krig_rmse(pred5, n)
pred6_rmse <- krig_rmse(pred6, n)
pred6.5_rmse <- krig_rmse(pred6.5, n)


#RMSE
rmse0 <- rmse(pred_rmse$elev_obs, pred_rmse$elev_pred) 
rmse2 <- rmse(pred2_rmse$elev_obs, pred2_rmse$elev_pred) 
rmse3 <- rmse(pred3_rmse$elev_obs, pred3_rmse$elev_pred) 
rmse3.5 <- rmse(pred3.5_rmse$elev_obs, pred3.5_rmse$elev_pred) 
rmse4 <- rmse(pred4_rmse$elev_obs, pred4_rmse$elev_pred) 
rmse5 <- rmse(pred5_rmse$elev_obs, pred5_rmse$elev_pred) 
rmse6 <- rmse(pred6_rmse$elev_obs, pred6_rmse$elev_pred) 
rmse6.5 <- rmse(pred6.5_rmse$elev_obs, pred6.5_rmse$elev_pred) 

#MSE
mse0 <- mse(pred_rmse$elev_obs, pred_rmse$elev_pred) 
mse2 <- mse(pred2_rmse$elev_obs, pred2_rmse$elev_pred) 
mse3 <- mse(pred3_rmse$elev_obs, pred3_rmse$elev_pred) 
mse3.5 <- mse(pred3.5_rmse$elev_obs, pred3.5_rmse$elev_pred) 
mse4 <- mse(pred4_rmse$elev_obs, pred4_rmse$elev_pred) 
mse5 <-mse(pred5_rmse$elev_obs, pred5_rmse$elev_pred) 
mse6 <- mse(pred6_rmse$elev_obs, pred6_rmse$elev_pred) 
mse6.5 <- mse(pred6.5_rmse$elev_obs, pred6.5_rmse$elev_pred) 


############################################################################################################
#CROSS VALIDATION
############################################################################################################

# Comparison against held out data (OBSERVED)
held_out_sample <- held_out   
held_out_sto <- stConstruct(held_out_sample, c("LAT", "LON"), "TIME", interval = FALSE)  
held_out_sto@sp@proj4string <- CRS('+init=epsg:3310')
held_out_sp <- as(held_out_sto, "Spatial")
time_list_obs <- unique(held_out_sp@data$time)  #full list of observations in time

krig_rmse <- function(krigged_data, validation_sample_size) {
  residuals <- list()
  predsp <- as(krigged_data, "Spatial") 
  gridded(predsp) <- TRUE
  time_list_pred <- unique(attr(predsp, "time"))
  
  for (i in 1:validation_sample_size) {
    
    #OBSERVED DATASET
    #select point in time
    random_time <- sample(1:length(time_list_obs), 1)  
    time_obs <- time_list_obs[random_time]   
    names(time_obs) <- "time_obs"
    #one observed point in space
    space_obs <- held_out_sp[held_out_sp$time == time_obs,]  
    random_space_obs <- sample(1:length(space_obs$ELEV), 1)
    space_obs <- space_obs[random_space_obs,]
    elev_obs <- space_obs$ELEV
    coor_obs <- space_obs@coords[random_space_obs,]
    lon <- coor_obs[1]
    names(lon) <- "lon"
    lat <- coor_obs[2]
    names(lat) <- "lat"
    names(elev_obs) <- "elev_obs"
    #PREDICTED DATASET
    #time
    time_pred = paste(substr(as.character(time_obs), 1,7),"-01", sep="")
    #space
    space_pred <- raster(predsp[time_pred])
    
    #EXTRACT PREDICTED RASTER DATA FROM OBSERVED POINT 
    elev_pred <- extract(space_pred, space_obs@coords)
    names(elev_pred) <- "elev_pred"
    
    #final data entry
    residuals[[i]] <- list(lon, lat, time_obs, elev_obs, time_pred, elev_pred)  #add time_obs, add grid.index coords for predicted
  }
  
  res <- as.data.frame(do.call(rbind, lapply(residuals, unlist)))
  res <- na.omit(res)
  res$time_pred <- as.Date(res$time_pred)
  res$time_obs <- as.Date(res$time_obs)
  
  return(res)  #return list, MSE, RMSE and residuals
}

#build near full res dataset
res <- krig_rmse(pred2, 47000)
# res$elev_predt<-log(res$elev_pred + constant_transform) #if reading in saved krigged and transformed dataset, back transform the predictions first
res$elev_delta <- res$elev_obs - res$elev_pred
hist(res$elev_delta, 1000)  

#Calucalte RMSE and Nash-Sutcliffe Efficiency
rmse<-rmse(res$elev_predt, res$elev_obs) #if reading in saved krigged & transformed data used backtransformed predictions
NSE<-NSE(res$elev_predt, res$elev_obs)
nrmse<-nrmse(res$elev_predt, res$elev_obs, norm = "maxmin")
meanNrmse<-rmse/mean(res$elev_obs, na.rm=TRUE)

#Build a space-time object for residuals data
#sto_res <- stConstruct(res, c("lat", "lon"), "time_pred", interval = FALSE)  #time instance, not time interval
#stplot(sto_res)

#plot residuals value at a point in space (held out) through all time obs
time_obs <- res%>% group_by(time_pred) %>% summarize(mean(elev_obs))
time_pred <- res%>% group_by(time_pred) %>% summarize(mean(elev_pred))
time_m <- res %>% group_by(time_pred) %>% summarize(mean(elev_delta))
time_sd <- res %>% group_by(time_pred) %>% summarize(sd(elev_delta))
time_count <- held_out_sample%>% group_by(TIME) %>% summarize(count=n())

ggplot(time_m, aes(time_m$time_pred, time_m$`mean(elev_delta)`)) +
  geom_line() + 
  labs(x = "Time", y="Observed - predicted elevation (feet)") +
  geom_ribbon(aes(ymin = time_m$`mean(elev_delta)`-time_sd$`sd(elev_delta)`,
                  ymax = time_m$`mean(elev_delta)` + time_sd$`sd(elev_delta)`), alpha=.2)


plot(time_m, type="l", main="Observed - predicted elevation", col="black", ylim=c(-10, 10), 
     ylab='Mean elevation difference (feet)', xlab='Time') 

#lines(time_obs, type='l', col='blue', lty=2)
#lines(time_pred, type='l', col='red', lty=2)
#ylab('Mean elevation difference (meters)')
#xlab('Time')
#lines(time_count, type='l', col='green')
#legend("bottomleft", c("obs", "pred", "difference"), lty=c(1,1), lwd=c(2.5,2.5), col=c("blue","red", "black")) 

#Plot residulas through space
space_m <- data.frame(res %>% group_by(lat,lon) %>% summarize(mean(elev_delta)))
space_sd <- data.frame(res %>% group_by(lat,lon) %>% summarize(sd(elev_delta)))
coordinates(space_m) <- ~lon+lat
space_m@proj4string <- CRS('+init=epsg:3310') #confirm original projection
coordinates(space_sd) <- ~lon+lat
space_sd@proj4string <- CRS('+init=epsg:3310') #confirm original projection
spplot(space_m, main="Observed - predicted elevation (feet)") #CHANGE SYMBOLOGY
cv_shp <- shapefile("C:\\Users\\Emily\\Box Sync\\GIS\\USA\\CA\\CA_state.shp")
cv_shp <- spTransform(cv_shp, space_m@proj4string)
plot(cv, add=T)
spplot(space_sd, main="obs-pred sd")

space_m_obs <- data.frame(res %>% group_by(lat,lon) %>% summarize(max(elev_obs)))
coordinates(space_m_obs) <- ~lon+lat
space_m_obs@proj4string <- CRS('+init=epsg:3310') #confirm original projection
spplot(space_m_obs, main="observed max")

space_m_p <- data.frame(res %>% group_by(lat,lon) %>% summarize(max(elev_pred)))
coordinates(space_m_p) <- ~lon+lat
space_m_p@proj4string <- CRS('+init=epsg:3310') #confirm original projection
spplot(space_m_p, main="predicted max")

#plot space-time fit
#extract all spatial observations at one point in time
select_time <- pred@data$var1.pred[pred@time$timeIndex == 6]

#plot elevation points through time
#predicted
grid_id <- unique(pred@sp@grid.index)
obs <- sample(grid_id, 1)
coords <- pred@sp@coords[pred@sp@grid.index == obs]
select_space <- pred@data$var1.pred[pred@sp@grid.index == obs]
plot(select_space, main=paste(obs, coords))

#observed
obs <- sample(unique(sto@sp@coords), 1)
coords <- sto@sp@coords[obs,]  
select_space <- sto@data$ELEV[sto@sp@coords == coords]
plot(select_space, main=paste(obs, coords))


#back transform krigged elevations
krigged_data <- pred
krigged_data@data <- exp(pred@data) - constant_transform
saveRDS(krigged_data, '/data/emily/WF/kate/final_data/krigged_data.rds')
kd <- readRDS('/data/emily/WF/kate/final_data/krigged_data.rds')
