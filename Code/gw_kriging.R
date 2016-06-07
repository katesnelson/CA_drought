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
#dr <- '/data/emily/WF/kate/gw/'
#d <- read.csv(paste(dr,'gw.csv',sep=''), stringsAsFactors=FALSE)

dr <- 'C:\\Users\\Emily\\Dropbox\\Vanderbilt\\Kate_Emily\\CA_drought\\Code\\'
d <- as.data.frame(read.table(paste(dr, 'dsub.txt', sep=''), header=T), stringsAsFactors=F)
coordinates(d) <- ~LON+LAT

#create space-time object, http://www.inside-r.org/packages/cran/spacetime/docs/stConstruct
df <- data.frame(log(d@data$ELEV), as.Date(d@data$DATE), as.data.frame(d@coords[,1]), as.data.frame(d@coords[,2]))
names(df) <- c("LOG_ELEV", "TIME", "LON", "LAT")
df$LOG_ELEV[is.na(df$LOG_ELEV)] <- 0
sto <- stConstruct(df, c("LAT", "LON"), "TIME", interval = FALSE)  #time instance, not time interval
sto@sp@proj4string <- CRS('+init=epsg:4326')  #WGS84
stplot(sto)

#space-time SAMPLE variogram
sampleVar <- variogramST(LOG_ELEV~1, sto, tunit="days", tlags=seq(from=0,to=225,by=30), cutoff = 100, na.omit=T)  #log elevation doesn't work because of negative values

#diagnostic plots
acf(sto@data, 'xts')  #up to 200 days significant
plot(sampleVar, map=F) #range of around 60, nugget less than 1
plot(sampleVar, wireframe=T) #visualize space/time lag
#for unprojected great circle distances (km).
#gamma values conditional on treatment of data (log, etc)

#estimate space-time anisotropy, or directional dependence of data
#http://finzi.psych.upenn.edu/library/gstat/html/estiStAni.html
estiStAni(sampleVar, c(10,150), "vgm", vgm(20,"Exp",120,0))
#data too sparse?  http://stats.stackexchange.com/questions/35316/problems-estimating-anisotropy-parameters-for-a-spatial-model

#parameter range
pars.l <- c(sill.s = 0, range.s = 10, nugget.s = 0, sill.t = 0, range.t = 10, 
            nugget.t = 0, sill.st = 0, range.st = 10, nugget.st = 0, anis = 0)
pars.u <- c(sill.s = 100, range.s = 100, nugget.s = 10,sill.t = 300, range.t = 100, 
            nugget.t = 100, sill.st = 200, range.st = 1000, nugget.st = 100,anis = 700) 

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
                joint = vgm(20,"Exp",100,0), #sill, model, range, nugget
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

#space-time-kriging document
#sampleVar
#plot empirical vars
#create version of table 2, figure 4 shows space-time biases
#find model with smallest RMSE

#plot all 2D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF), all=T) #psVarBF add once k issue fixed

#plot all 3D
plot(sampleVar, list(metVarBF, smVarBF, ssmVarBF), all=T, wireframe=T, zlim=c(0,120),  #psVarBF add once k issue fixed
     zlab=NULL, xlab=list("distance (km)", rot=30),
     ylab=list("time lag (days)", rot=-35),
     scales=list(arrows=F, z = list(distance = 5)))

#plot fit metrics (data fit)
MSE_DF <- c(MSEmetVarDF, MSEsmVarDF, MSEssmVarDF) #add ps
barplot(MSE_DF, main="MSE of data fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Simple-sum-metric"))

#plot fit metrics (best fit)
MSE_BF <- c(MSEmetVarBF, MSEsmVarBF, MSEssmVarBF)  #add ps
barplot(MSE_BF, main="MSE of best fit", ylab="MSE", xlab = "Variogram model",
        names.arg=c("Metric", "Simple-metric", "Simple-sum-metric"))


#create new spatiotemporal prediction grid
pred<-krigeST(LOG_ELEV~1, data=sto, newdata=NEWDATAGRID, sepVar, nmax=50)
stplot(pred)


#leave-one-out cross-validation
#plot figure 6, differences at locations, and 7, differences through time




#fit tables:
MSEmetVarDF <- attr(metVarDF, "MSE")
MSEmetVarBF <- attr(metVarBF, "MSE") 
