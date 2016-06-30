#questions/issues:
#should we do any additional pre-processing of the elevation data (standardize, log, etc)?
#space-time SAMPLE variogram, does not work with neg values

#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/
#http://localhost/rstudio//files/R/x86_64-pc-linux-gnu-library/3.2/gstat/doc/st.pdf

library(dplyr)
library(gstat)
library(spacetime) #https://cran.r-project.org/web/packages/spacetime/vignettes/jss816.pdf
library(sp)  #creates spatial object
library(xts) #creates time object
library(reshape2)
library(ggplot2)
library(rgdal)
library(raster)

#ekb comp
#dr <- 'C:\\Users\\Emily\\Dropbox\\Vanderbilt\\Kate_Emily\\CA_drought\\Code\\'
#d <- as.data.frame(read.table(paste(dr, 'dsub.txt', sep=''), header=T), stringsAsFactors=F)

#ornette
d <- as.data.frame(read.csv(paste('/data/emily/WF/kate/gw/gw.csv', sep=','), header=T), stringsAsFactors=F)
coordinates(d) <- ~LON+LAT

#load data
df <- data.frame(d@data$ELEV, as.Date(d@data$DATE), as.data.frame(d@coords[,1]), as.data.frame(d@coords[,2]))
names(df) <- c("ELEV", "TIME", "LON", "LAT")

#clean up outliers, drop values less than -50 (10%) and more than 3000 (2%)
df <- subset(df, df$ELEV >= -50 & df$ELEV < 3000, select=c("ELEV", "TIME", "LON", "LAT"))

#select random subset of data (10%)
dfs <- df[sample(nrow(df), 0.1*nrow(df)), ]

#space-time object
sto <- stConstruct(dfs, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval
sto@sp@proj4string <- CRS('+init=epsg:4326')  #WGS84
stplot(sto)

#create sample variogram, NOTE: takes awhile to run
#sampleVar <- variogramST(ELEV~1, sto, tunit="days", tlags=seq(from=0,to=225,by=15), cutoff = 40, na.omit=T)  #log elevation doesn't work because of negative values
#write.table(sampleVar, 'sampleVar.csv', sep=",")
sampleVar <- read.table('sampleVar.csv', sep=",")

  #avg width of Central Valley is 64 km, 0.58*64 is 37; round up to 40 for spatial range
  #acf(sto@data, 'xts')  #up to 200 days significant
  #acf is nonstationary, but most change occurs quickly within the first month, so smaller interval should be used, too small and you get holes
  #from acf most significance within 6 months or 180 days 


#diagnostic plots
plot(sampleVar, map=F) 
plot(sampleVar, wireframe=T) 
#NOTE:  gamma values conditional on treatment of data (log, etc)


#parameter range, http://giv-graeler.uni-muenster.de:3838/spacetime/; provide huge range here
pars.l <- c(sill.s = 0.1, range.s = 10, nugget.s = 0.1, sill.t = 0.1, range.t = 10, 
            nugget.t = 0.1, sill.st = 0.1, range.st = 10, nugget.st = 0.1, anis = 0.1)
pars.u <- c(sill.s = 500000, range.s = 200000, nugget.s = 3000, sill.t = 300000, range.t = 200000, 
            nugget.t = 1000, sill.st = 3000, range.st = 100000, nugget.st = 10000, anis = 70000) 

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
                joint = vgm(100,"Exp",200,0), #sill, model, range, nugget
                stAni = 100) #directional dependence
metVarDF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=0)
MSEmetVarDF <- attr(metVarDF, "MSE")
metVarBF <- fit.StVariogram(sampleVar, metVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEmetVarBF <- attr(metVarBF, "MSE") 

#fit.model = 7 is method that controls the weighing of the residuals between empirical and model surface;
#7 is the ST analog to commonly used spatial weighting; higher confidence in distances with many pairs

###################
#product-sum
##################
#relax assummptions in covariance matrix, does not assume separability
#need to provude both spatial and temporal components, plus k, weighting product

psVar <- vgmST("productSum", 
                space=vgm(15, "Exp", 60, 0),   #sill, model, range, nugget
                time=vgm(15, "Exp", 40, 0),
                k=100)
psVarDF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=0)  #data fit
MSEpsVarDF <- attr(psVarDF, "MSE")
psVarBF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
#k must be positive error
MSEpsVarBF <- attr(psVarBF, "MSE")

####################
#sum-metric
####################
#metric with spatial and temporal covariance models plus joint component with anisitropy
#maximum flexibility

smVar <- vgmST("sumMetric", 
                space=vgm(15, "Exp", 60, 0),   #sill, model, range, nugget
                time=vgm(15, "Exp", 40, 0),
                joint=vgm(15,"Exp",80,0),
                stAni=200)
smVarDF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=0)  
MSEsmVarDF <- attr(smVarDF, "MSE")
smVarBF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEsmVarBF <- attr(smVarBF, "MSE")

######################
#simple sum-metric
######################
#like sum metric but restrict all components to have a single nugget

ssmVar <- vgmST("simpleSumMetric", 
                space=vgm(15, "Exp", 60, 0),   #sill, model, range, nugget
                time=vgm(15, "Exp", 40, 0),
                joint=vgm(15,"Exp",80,0),
                nugget=1, stAni=77)
ssmVarDF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=0)  
MSEssmVarDF <- attr(ssmVarDF, "MSE")
ssmVarBF <- fit.StVariogram(sampleVar, ssmVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
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

#product sum wins

#############################################################
# Prediction grid construction
#############################################################

#create empty spatial object
so <- raster()
extent(so) <- extent(d)
#8.9858 degrees x by 7.13729 degrees y, 1 degree, 110 km... approx 2.75 km cell size
so <- as(so, 'SpatialPoints') #360 x 180 cells
gridded(so) <- TRUE
so@proj4string <- sto@sp@proj4string  #set sto proj

#create monthly time object
tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 12*15)

#sto frame
grid.st <- STF(sp = as(so, "SpatialPoints"), time = tm.grid)

#create new spatiotemporal prediction grid
pred <- krigeST(ELEV~1, data=sto, newdata=grid.st, modelList=psVarBF, bufferNmax=1, progress=T)
gridded(pred@sp) <- TRUE 
stplot(pred)

#leave-one-out cross-validation
#plot figure 6, differences at locations, and 7, differences through time

