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

#jg direction in space-time variography
#what's up with variogram bumps
#any preprocessing of groundwater elevtaion data?  detrend?  log

#home comp
#dr <- 'C:\\Users\\Emily\\Dropbox\\Vanderbilt\\Kate_Emily\\CA_drought\\Code\\'
#d <- as.data.frame(read.table(paste(dr, 'dsub.txt', sep=''), header=T), stringsAsFactors=F)

#ornette
d <- as.data.frame(read.csv(paste('/data/emily/WF/kate/gw/gw.csv', sep=','), header=T), stringsAsFactors=F)
coordinates(d) <- ~LON+LAT

#create space-time object, http://www.inside-r.org/packages/cran/spacetime/docs/stConstruct
#take log of elev for more reasonable numbers, but negative and near-zero numbers cause NANs
df <- data.frame(d@data$ELEV, as.Date(d@data$DATE), as.data.frame(d@coords[,1]), as.data.frame(d@coords[,2]))
names(df) <- c("ELEV", "TIME", "LON", "LAT")

#clean up outliers, drop values less than -50 (10%) and more than 3000 (2%)
df <- subset(df, df$ELEV >= -50 & df$ELEV < 3000, select=c("ELEV", "TIME", "LON", "LAT"))

#random subset of data
dfs <- df[sample(nrow(df), 0.02*nrow(df)), ]
dfs<- distinct(dfs, TIME, LAT, LON, .keep_all=TRUE)

#space-time object
sto <- stConstruct(dfs, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval
sto@sp@proj4string <- CRS('+init=epsg:4326')  #WGS84
stplot(sto)

#space-time SAMPLE variogram, does not work with neg values
#avg width of CA Central Valley is 64 km, 0.58*64 is 37 --> round up to 40 for spatial range
#acf is nonstationary, but most change occurs quickly within the first month, so smaller interval should be used, too small and you get holes
#from acf most significance within 6 months or 180 days 
sampleVar <- variogramST(ELEV~1, sto, tunit="days", tlags=seq(from=0,to=225,by=15), cutoff = 60, width = 5,  na.omit=T)  #trying to specify a small width to see if we can resolve possible high correlation issues
#write.table(sampleVar, 'sampleVar.csv', sep=",")
#sampleVar <- read.table('sampleVar.csv', sep=",")

#diagnostic plots
acf(sto@data, 'xts')  #up to 200 days significant
plot(sampleVar, map=F) #range of around 60, nugget less than 1
plot(sampleVar, wireframe=T) #visualize space/time lag
#for unprojected great circle distances (km).
#gamma values conditional on treatment of data (log, etc)

#estimate space-time anisotropy, or directional dependence of data
#http://finzi.psych.upenn.edu/library/gstat/html/estiStAni.html
t<-estiStAni(sampleVar, interval = c(1,50), method = "metric", vgm(550000, "Gau", 46, 43999), s.range = 75, t.range = 180)
#data too sparse?  http://stats.stackexchange.com/questions/35316/problems-estimating-anisotropy-parameters-for-a-spatial-model

#parameter range
#paramcheck:http://giv-graeler.uni-muenster.de:3838/spacetime/
pars.l <- c(sill.s = 0.1, range.s = 10, nugget.s = 0.1, sill.t = 0.1, range.t = 1, 
            nugget.t = 0.1, sill.st = 0.1, range.st = 10, nugget.st = 0.1, anis = 0.0001)
pars.u <- c(sill.s = 2000000, range.s = 150, nugget.s = 200000, sill.t = 500000, range.t = 365, 
            nugget.t = 200000, sill.st = 3000000, range.st = 100000, nugget.st = 100000, anis = 3) #upper anis based on estiStAni

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
                joint = vgm(350000,"Gau",150,180000), #sill, model, range, nugget
                stAni = 500000) #directional dependence
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
               space=vgm(25000, "Gau", 42, 44),   #sill, model, range, nugget
               time=vgm(185, "Exp", 365, 100),
               k=9000)
psVarDF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=0)  #data fit
MSEpsVarDF <- attr(psVarDF, "MSE")
psVarBF <- fit.StVariogram(sampleVar, psVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u, stAni = .4)
#k must be positive error
MSEpsVarBF <- attr(psVarBF, "MSE")

####################
#sum-metric
####################
#metric with spatial and temporal covariance models plus joint component with anisitropy
#maximum flexibility

smVar <- vgmST("sumMetric", 
               space=vgm(550000, "Gau", 45, 44000),   #sill, model, range, nugget
               time=vgm(317000, "Exp", 10, 37000),
               joint=vgm(37000,"Gau",80, 46000),
               stAni=70000)
smVarDF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=0)  
MSEsmVarDF <- attr(smVarDF, "MSE")
smVarBF <- fit.StVariogram(sampleVar, smVar, method="L-BFGS-B", fit.method=7, lower=pars.l, upper=pars.u)
MSEsmVarBF <- attr(smVarBF, "MSE")

######################
#simple sum-metric
######################
#like sum metric but restrict all components to have a single nugget

ssmVar <- vgmST("simpleSumMetric", 
                space=vgm(750000, "Gau", 60),   #sill, model, range, nugget
                time=vgm(0.1, "Exp", 0.1),
                joint=vgm(11,"Gau",80,0),
                nugget=350000, stAni=100000)
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

#plot fit metrics (data fit)
MSE_DF <- c(MSEmetVarDF, MSEsmVarDF, MSEssmVarDF, MSEpsVarDF) 
barplot(MSE_DF, main="MSE of data fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Sum-metric", "Product-sum"))

#plot fit metrics (best fit)
MSE_BF <- c(MSEmetVarBF, MSEsmVarBF, MSEssmVarBF, MSEpsVarBF)  
barplot(MSE_BF, main="MSE of best fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Sum-metric", "Product-sum"))

#extractPar(smVarBF) #gives fitted parameter values, anis is always the max?
#############################################################
# Prediction grid construction
#############################################################

#create spatial reference
newproj <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
rast <- raster()
extent(rast) <- extent(d)
rast <- as(rast, 'SpatialPixels') #360 x 180 cells
res(rast@grid@cells.dim) <- c(180,90)


#create time object
#tm.grid <- dfs$TIME[order(dfs$TIME)]
tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 12*15)

tm.grid <- seq(as.Date("2000/1/1"), by = "month", length.out = 3)
#sto frame
grid.st <- STF(rast, tm.grid)

sto@sp@proj4string <- CRS('+init=epsg:3310')  # set projection for kriging to NAD83/california albers
grid.st=spTransform(grid.st,sto@sp@proj4string) #set to same proj as sto


#create new spatiotemporal prediction grid
pred<-krigeST(ELEV~1, data=sto, newdata=grid.st, nmax = 20, modelList=smVarBF, stAni = 1, progress=TRUE) #specifiying an nmax reduces runtime, but requires and stAni
stplot(pred)

#http://stats.stackexchange.com/questions/110156/r-error-in-chol-defaulta-krigest-from-gstat-package  
# A nugget of 0 may lead to Error in chol.default(A) : the leading minor of order 17 is not positive definite
# Check of variogram model indicates that the nugget is not 0

#http://stackoverflow.com/questions/25063157/r-error-in-chol-defaulta-krigest-from-gstat-package  
# duplicate observations or variogram that does not discrimiate observations sufficiently, leading to near perfect correlation of obs
# We don't appear to have duplicate observations (at least in the subset I pulled)
# Problem is probably near perfect correlation between obs, 
# maybe can be fixed by reducing the width of the variogram model spatial increments or by using a smaller set of obs



#leave-one-out cross-validation
#plot figure 6, differences at locations, and 7, differences through time

