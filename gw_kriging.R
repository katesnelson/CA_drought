#http://www.r-bloggers.com/spatio-temporal-kriging-in-r/
#http://localhost/rstudio//files/R/x86_64-pc-linux-gnu-library/3.2/gstat/doc/st.pdf

library(dplyr)
library(gstat)
library(spacetime)
library(sp)
library(reshape2)

dr <- '/data/emily/WF/kate/gw/'

d <- read.csv(paste(dr,'gw.csv',sep=''), stringsAsFactors=FALSE)

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