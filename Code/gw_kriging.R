#to do: rerun sampleVar and save

#questions/issues:
#should we do any additional pre-processing of the elevation data (detrend, standardize, log, etc)?
#space-time SAMPLE variogram, does not work with neg values
#check initial projection definition with kate's data
#ok to drop negative log values?

#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/
#http://localhost/rstudio//files/R/x86_64-pc-linux-gnu-library/3.2/gstat/doc/st.pdf

library(dplyr)
library(gstat)
library(spacetime) 
library(sp)  
library(xts) 
library(rgdal)
library(raster)
library(Metrics)

set.seed(10010)

#load data on ornette
d <- as.data.frame(read.csv(paste('/data/emily/WF/kate/gw/gw.csv', sep=','), header=T), stringsAsFactors=F)
coordinates(d) <- ~LON+LAT
df <- data.frame(d@data$ELEV, as.Date(d@data$DATE), as.data.frame(d@coords[,1]), as.data.frame(d@coords[,2]))
names(df) <- c("ELEV", "TIME", "LON", "LAT")

#clean up outliers, drop values less than 0 (10%) and more than 2500 (10%)
df <- subset(df, df$ELEV >= 0.001 & df$ELEV < 2500, select=c("ELEV", "TIME", "LON", "LAT"))
df$ELEV <- log(df$ELEV)

#select random subset of data (10%)
#does not necessarily cover all space-time points with observations
dfs <- df[sample(nrow(df), 0.1*nrow(df)), ] 

#remove any duplicates
dfs_clean <- dfs[!(duplicated(dfs[c("LAT", "LON", "TIME")]) | duplicated(dfs[c("LAT","LON", "TIME")], fromLast = TRUE)), ]

#sample of held out datapoints to test fit later
held_out <- rbind(df, dfs) #duplicate sampled observations
held_out <- held_out[! duplicated(held_out, fromLast=TRUE) & seq(nrow(held_out)) <= nrow(df), ]

#clean up unnecessary files
remove(d, df, dfs)

#space-time object
sto <- stConstruct(dfs_clean, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval
sto@sp@proj4string <- CRS('+init=epsg:3310')  # set projection for kriging to NAD83/california albers
#stplot(sto)

#create sample variogram
sampleVar <- variogramST(ELEV~1, sto, tunit="days", tlags=seq(from=0,to=180,by=15), cutoff = 60, na.omit=T)  
#write.table(sampleVar, 'sampleVar.csv', sep=",")
#sampleVar <- read.table('sampleVar.csv', sep=",")

#diagnostic plots
plot(sampleVar, map=F) 
plot(sampleVar, wireframe=T) 

#parameter range, http://giv-graeler.uni-muenster.de:3838/spacetime/
#provide huge range here
pars.l <- c(sill.s = 0.1, range.s = 10, nugget.s = 0.1, sill.t = 0.1, range.t = 10, 
            nugget.t = 0.1, sill.st = 0.1, range.st = 10, nugget.st = 0.1, anis = 0.1)
pars.u <- c(sill.s = 500000, range.s = 200000, nugget.s = 3000, sill.t = 300000, range.t = 200000, 
            nugget.t = 1000, sill.st = 3000, range.st = 100000, nugget.st = 10000, anis = 70000) 



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

metVar <- vgmST("metric",
                joint = vgm(100,"Exp",200,0), #sill, model, range, nugget
                stAni = 100) #directional dependence
metVarDF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=0) #data fit
MSEmetVarDF <- attr(metVarDF, "MSE")
metVarBF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEmetVarBF <- attr(metVarBF, "MSE") 

############################################################################################################
#PRODUCT SUM
############################################################################################################

#relax assummptions in covariance matrix, does not assume separability
#need to provide both spatial and temporal components, plus k, weighting product

psVar <- vgmST("productSum", 
                space=vgm(15, "Exp", 60, 0), 
                time=vgm(15, "Exp", 40, 0),
                k=100)
psVarDF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=0)  
MSEpsVarDF <- attr(psVarDF, "MSE")
psVarBF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEpsVarBF <- attr(psVarBF, "MSE")


############################################################################################################
#SUM METRIC
############################################################################################################

#metric with spatial and temporal covariance models plus joint component with anisitropy
#maximum flexibility

smVar <- vgmST("sumMetric", 
                space=vgm(15, "Exp", 60, 0),
                time=vgm(15, "Exp", 40, 0),
                joint=vgm(15,"Exp",80,0),
                stAni=200)
smVarDF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=0)  
MSEsmVarDF <- attr(smVarDF, "MSE")
smVarBF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEsmVarBF <- attr(smVarBF, "MSE")

############################################################################################################
#SIMPLE SUM METRIC
############################################################################################################

#like sum metric but restrict all components to have a single nugget

ssmVar <- vgmST("simpleSumMetric", 
                space=vgm(15, "Exp", 60, 0),
                time=vgm(15, "Exp", 40, 0),
                joint=vgm(15,"Exp",80,0),
                nugget=1, stAni=77)
ssmVarDF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=0)  
MSEssmVarDF <- attr(ssmVarDF, "MSE")
ssmVarBF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEssmVarBF <- attr(ssmVarBF, "MSE")

############################################################################################################
#SEPARABLE MODEL
############################################################################################################

#not using since bad assumptions, but documentation and scripts below:
#http://geostat-course.org/system/files/part01.pdf: assumptions of isotropy and stationarity
#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/: 
#pararmeter selection for vgmST doesn't seem to matter so much 

# sepVar <- vgmST("separable", 
#                   space=vgm(15, "Exp", 60, 0),   #sill, model, range, nugget
#                   time=vgm(15, "Exp", 40, 0),
#                   sill = 30,
#                   stAni = 100)
# sepVarDF <- fit.StVariogram(sampleVar, sepVar, method="L-BFGS-B", fit.method=0)  
# MSEsepVarDF <- attr(sepVarDF, "MSE")  
# sepVarBF <- fit.StVariogram(sampleVar, sepVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
# MSEsepVarBF <- attr(sepVarBF, "MSE") 






############################################################################################################
#MODEL COMPARISON
############################################################################################################

#plot all 2D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF, psVarBF), all=T) 

#plot all 3D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF, psVarBF), all=T, wireframe=T, 
     zlab=NULL, xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)))

#plot fit metrics (data fit)
MSE_DF <- c(MSEmetVarDF, MSEsmVarDF, MSEssmVarDF, MSEpsVarDF) 
barplot(MSE_DF, main="MSE of data fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Sum-metric", "Product-sum"))

#plot fit metrics (best fit)
MSE_BF <- c(MSEmetVarBF, MSEsmVarBF, MSEssmVarBF, MSEpsVarBF)  
barplot(MSE_BF, main="MSE of best fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Sum-metric", "Product-sum"))





############################################################################################################
#PREDICTION GRID CONSTRUCTION
############################################################################################################

#create empty spatial object
so <- raster()
extent(so) <- extent(d)
so <- setValues(so, 0)
so_coarse <- aggregate(so, fact=10, fun=mean, expand=TRUE, na.rm=TRUE)
so_coarse@crs <- sto@sp@proj4string  #set to sto proj

#create monthly time object
tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 12*15)

#empty sto object
grid.st <- STF(sp = as(so_coarse, "SpatialPixels"), time = tm.grid)

#krig data
pred <- krigeST(ELEV~1, data=sto, newdata=grid.st, nmax=20, modelList=psVarBF, stAni=1, progress=T)
gridded(pred@sp) <- TRUE 
stplot(pred)

#add stAni, change nmax, manipulate grid, change model, etc.





############################################################################################################
#CROSS VALIDATION
############################################################################################################

#random sample of held-out data (OBSERVED)
held_out_sample <- held_out[sample(nrow(held_out), 0.01*nrow(held_out)), ]
held_out_sto <- stConstruct(held_out_sample, c("LAT", "LON"), "TIME", interval = FALSE)  
held_out_sto@sp@proj4string <- CRS('+init=epsg:3310')
held_out_sp <- as(held_out_sto, "Spatial")
time_list_obs <- unique(held_out_sp@data$time)

#rasterize predicted sto object (PREDICTED)
predsp <- as(pred, "Spatial") #spatial pixels dataframe
gridded(predsp) <- TRUE
time_list_pred <- unique(attr(predsp, "time"))


validation_sample_size <- 100
residuals <- list()

for (i in 1:validation_sample_size) {
  
  #predicted data
  random_time <- sample(1:length(time_list_pred), 1)
  time_pred <- time_list_pred[random_time]
  names(time_pred) <- "time_pred"
  space_pred <- raster(predsp[random_time])

  #observed data
  time_gt <- time_list_obs[time_list_obs > time_pred]
  if (length(time_gt) == 0) {
    time_lt <- time_list_obs[time_list_obs < time_pred]
    time_obs <- tail(time_lt, n=1)
  } else {
    time_obs <- time_gt[1]
  }
  space_obs <- held_out_sp[held_out_sp$time == time_obs,]
  elev_obs <- space_obs$ELEV[1]
  coor <- space_obs@coords[1,]
  lon <- coor[1]
  names(lon) <- "lon"
  lat <- coor[2]
  names(lat) <- "lat"
  names(elev_obs) <- "elev_obs"
  
  #extract predicted data from observed data point
  elev_pred <- extract(space_pred, space_obs)[1]
  names(elev_pred) <- "elev_pred"
  
  #final data entry
  residuals[[i]] <- list(elev_obs, elev_pred, time_pred, lon, lat)
}

res <- as.data.frame(do.call(rbind, lapply(residuals, unlist)))
res <- na.omit(res)
res$time_pred <- as.Date(res$time_pred)

mse <- mse(res$elev_obs, res$elev_pred)
rmse <- rmse(res$elev_obs, res$elev_pred)
