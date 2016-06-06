#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/
#http://localhost/rstudio//files/R/x86_64-pc-linux-gnu-library/3.2/gstat/doc/st.pdf

library(dplyr)
library(gstat)
library(spacetime) #https://cran.r-project.org/web/packages/spacetime/vignettes/jss816.pdf
library(sp)  #creates spatial object
library(xts) #creates time object
library(reshape2)

dr <- '/data/emily/WF/kate/gw/'

d <- read.csv(paste(dr,'gw.csv',sep=''), stringsAsFactors=FALSE)
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