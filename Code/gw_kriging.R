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

#jg direction in space-time variography
#what's up with variogram bumps
#any preprocessing of groundwater elevtaion data?  detrend?  log

#random sample of full dataset
#build sample vg
#id best modeled vg

#dr <- '/data/emily/WF/kate/gw/'
#d <- read.csv(paste(dr,'gw.csv',sep=''), stringsAsFactors=FALSE)

dr <- 'C:\\Users\\Emily\\Dropbox\\Vanderbilt\\Kate_Emily\\Code\\'
d <- as.data.frame(read.table(paste(dr, 'dsub.txt', sep=''), header=T), stringsAsFactors=F)
coordinates(d) <- ~LON+LAT

#create space-time object, http://www.inside-r.org/packages/cran/spacetime/docs/stConstruct
df <- data.frame(log(d@data$ELEV), as.Date(d@data$DATE), as.data.frame(d@coords[,1]), as.data.frame(d@coords[,2]))
names(df) <- c("LOG_ELEV", "TIME", "LON", "LAT")
df$LOG_ELEV[is.na(df$LOG_ELEV)] <- 0
sto <- stConstruct(df, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval
sto@sp@proj4string <- CRS('+init=epsg:4326')  #WGS84
stplot(sto)

#acf
acf(sto@data, 'xts')  #up to 200 days significant

#space-time SAMPLE variogram
sampleVar <- variogramST(LOG_ELEV~1, sto, tunit="days", tlags=seq(from=0,to=225,by=30), cutoff = 100, na.omit=T)  #log elevation doesn't work because of negative values
plot(sampleVar, map=F)
plot(sampleVar, wireframe=T)
#for unprojected great circle distances (km).
#gamma values conditional on treatment of data (log, etc)

#parameters
pars.l <- c(sill.s = 0, range.s = 10, nugget.s = 0,sill.t = 0, range.t = 1, 
            nugget.t = 0,sill.st = 0, range.st = 10, nugget.st = 0, anis = 0)
pars.u <- c(sill.s = 100, range.s = 1000, nugget.s = 100,sill.t = 200, range.t = 60, 
            nugget.t = 100,sill.st = 200, range.st = 1000, nugget.st = 100,anis = 700) 
modelVar <- vgmST("sumMetric", 
                      space=vgm(20, "Exp", 150, 0.1),   #sill, model, range, nugget
                      time=vgm(20, "Exp", 100, 0.1), 
                      joint=vgm(20,"Exp",150,0),
                      stAni=15)

#check parameters
extractPar(modelVar)
plot(sampleVar, modelVar, map=F)

#test fit
comp <- fit.StVariogram(sampleVar, modelVar, fit.method=0)
attr(comp, "MSE")




fit.StVariogram(sampleVar, vgm(c("Exp", "sph")))

#variogram MODELS



#metric model with st anisotropy
metricVarJoint <- vgmST("metric", 
                        joint=vgm(psill=10,"Mat",range=100, nugget=1), stAni=50)

metricVarSep <- vgmST("separable", 
                      space=vgm(-60, "Sph", 500, 1), 
                      time=vgm(35, "Sph", 500, 1), sill = 100)

#vgm(psill, model, range, nugget, add.to, annis, kappa)

metricVar <- fit.StVariogram(sampleVar, metricVarSep)
#RMS-difference between surfaces of metric and sampleVar
attr(metricVar,"optim")$value 


var <- variogramST(ELEV~LAT+LON, data=sto, tunit="days", tlags=seq(from=0, to=225, by=15), assumeRegular=F, na.omit=T)
plot(var, map=F)
#assuming separability, simple model

#metric vgm
metricVgm <- vgmST("metric", joint=vgm(50,"Exp", 100,0), stAni=50)
metricVgm <- fit.StVariogram(var, metricVgm)
plot(var, metricVgm)
#


#lower and upper limits for variogram parameters, now from example directly
pars.l <- c(sill.s = 0, range.s = 10, nugget.s = 0,sill.t = 0, range.t = 1, 
            nugget.t = 0,sill.st = 0, range.st = 10, nugget.st = 0, anis = 0)
pars.u <- c(sill.s = 200, range.s = 1000, nugget.s = 100,sill.t = 200, range.t = 60, 
            nugget.t = 100,sill.st = 200, range.st = 1000, nugget.st = 100,anis = 700)


#for separable variogram model, create model for spatial AND temporal components
sep <- vgmST("separable", space=vgm(-60, "Sph", 500, 1), time=vgm(35, "Sph", 500, 1), sill = 100)
par(mfrow=c(1,2))


#variogram plot
par(mfrow=c(1,2))
plot(var, map=F)
plot(var, sep, map=F)













#long format: each record reflects single space and time combination

#create space-time object, http://www.inside-r.org/packages/cran/spacetime/docs/stConstruct
d@data$LON <- d@coords[,1]
d@data$LAT <- d@coords[,2]
d@data$TIME <- as.Date(d@data$DATE)
sto <- stConstruct(d@data, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval

#variogram, http://www.inside-r.org/node/209774
var <- variogramST(ELEV~1, locations=data@coords, data=sto, tunit="days", tlags=seq(from=15, to=90, by=15), assumeRegular=F, na.omit=T, progress=T)

plot(var, map=F) #gama distance
plot(var, map=T) #time distance
plot(var, wireframe=T) #3d gamma distance time

#http://geostat-course.org/system/files/part01.pdf
tmp_pred <- data.frame(cbind(d@coords, d@data$TIME))
colnames(tmp_pred) <- c("x", "y", "t")
coordinates(tmp_pred) <- ~x+y+t
blockKrige <- krige(ELEV~1, sto, newdata=tmp_pred, model=model3d) #add block?



#variogram assumptions
Instead of using a different variogram for each day, we assume a constant model and use one variogram relying on
all data treating each day as a copy of the same process (pooled variogram) or averaging the variograms (mean
                                                variogram). The one estimated variogram is then used for
each day. This way, all information is used but again no
temporal interactions are incorporated.





#1.  Create SpatialPointsDataFrame

coordinates(d) =~ LON+LAT
projection(d)=CRS('+init=epsg:4326')  #WGS84
#transform to Mercator to read variogram in meters not degrees
d <- spTransform(d,CRS('+init=epsg:3395'))  

#2.  Create STIDF object (unstructured space-time objects)

dsp <- SpatialPoints(d@coords,CRS('+init=epsg:3395'))

#check for duplicate points
dupl <- zerodist(dsp)
df <- data.frame(d[-dupl[,2]])





tm <- as.POSIXct(d$DATE[-dupl[,2]])
timeDF <- STIDF(dsp,tm,data=df)

var <- variogramST(ELEV~1, data=timeDF, tunit="months", assumeRegular=F, na.omit=T)





#182 observations, 2000.01 to 2014.12, monthly

quality <- apply(data, MARGIN=1, FUN = function(x) length(x[!is.na(x)]))
#drop wells with fewer than 14 obs?
#n = 14345 wells total
#n = 5989 with > 14 observations

#lat/lon