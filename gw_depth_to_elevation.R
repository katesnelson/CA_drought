library(raster)
library(elevatr)


#This script extracts the groundwater elevation in January from the kriged prediction grid and returns a depth to groundwater by comparison to ground elevations. 

cv_huc <- shapefile('/data/emily/WF/kate/shp/cv_huc.shp')
#load krigged sto
kd <- readRDS('/data/emily/WF/kate/final_data/krigged_data.rds')

#from feet to meters
kd@data <- kd@data/3.28

#create elevation data
ras <- raster(kd[,1])
ext <- as(extent(ras), 'SpatialPolygons')  
ext@proj4string <- crs(ras)
elev <- get_elev_raster(ext, z=3, src="aws")  #zoom selected based on latitude and resolution of krigged dataset which is 10K
elev <- crop(elev, ext)

elev_rs <- resample(elev, ras, method="bilinear")

jan <- kd[ ,seq(1,length(kd@time),12)]
jan15 <- kd[,180]

elev_rs <- crop(elev_rs, raster(jan[,1]))

y01 <- elev_rs - raster(jan[,1])
y02 <- elev_rs - raster(jan[,2])
y03 <- elev_rs - raster(jan[,3])
y04 <- elev_rs - raster(jan[,4])
y05 <- elev_rs - raster(jan[,5]) #2K
y06 <- elev_rs - raster(jan[,6]) #-2K
y07 <- elev_rs - raster(jan[,7])
y08 <- elev_rs - raster(jan[,8])
y09 <- elev_rs - raster(jan[,9]) #2K
y10 <- elev_rs - raster(jan[,10])
y11 <- elev_rs - raster(jan[,11]) #1K
y12 <- elev_rs - raster(jan[,12])
y13 <- elev_rs - raster(jan[,13])
y14 <- elev_rs - raster(jan[,14])

diff <- brick(y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14)
cv_huc <- spTransform(cv_huc, diff@crs)
diff_crop <- crop(diff, cv_huc)
diff_val <- getValues(diff_crop)
hist(diff_val, 1000)

#drop values below zero
diff[diff < 0] <- NA
diff[diff > 2500] <- NA
names(diff) <- c("gw_01", "gw_02", "gw_03", "gw_04", "gw_05", "gw_06", "gw_07", "gw_08", "gw_09", "gw_10", "gw_11",
                 "gw_12", "gw_13", "gw_14")


d <- as.data.frame(read.csv('/data/emily/WF/kate/final_data/pt_data_May_05.txt', sep=','), header=T, stringsAsFactors=F)
coordinates (d)<- c("lon", "lat")
columns<-names(d)
d@proj4string<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

#join difference brick (diff) and point dataset (d)
d_prj <- spTransform(d, crs(diff)) #overlays on brick

#extract values of diff brick to points
diff_pt <- extract(diff, d_prj, method="simple", sp=T)

#formatting final dataset, from factor to numeric
diff_pt@data$POINTID <- as.numeric(diff_pt@data$POINTID)
diff_pt@data$FID <- as.numeric(diff_pt@data$FID)
diff_pt@data$GRID_CODE <- as.numeric(diff_pt@data$GRID_CODE)
diff_pt@data$RASTERVALU <- as.numeric(diff_pt@data$RASTERVALU)
diff_pt@data$spi_y3 <- as.numeric(as.character(diff_pt@data$spi_y3))  
diff_pt@data$spi_y4 <- as.numeric(as.character(diff_pt@data$spi_y4))
diff_pt@data$spi_y5 <- as.numeric(as.character(diff_pt@data$spi_y5))
diff_pt@data$spi_y6 <- as.numeric(as.character(diff_pt@data$spi_y6))
diff_pt@data$spi_y7 <- as.numeric(as.character(diff_pt@data$spi_y7))
diff_pt@data$spi_y8 <- as.numeric(as.character(diff_pt@data$spi_y8))
diff_pt@data$spi_y9 <- as.numeric(as.character(diff_pt@data$spi_y9))
diff_pt@data$spi_y10 <- as.numeric(as.character(diff_pt@data$spi_y10))
diff_pt@data$spi_y11 <- as.numeric(as.character(diff_pt@data$spi_y11))
diff_pt@data$spi_y12 <- as.numeric(as.character(diff_pt@data$spi_y12))
diff_pt@data$spi_y13 <- as.numeric(as.character(diff_pt@data$spi_y13))
diff_pt@data$spi_y14 <- as.numeric(as.character(diff_pt@data$spi_y14))
diff_pt@data$spi_y15 <- as.numeric(as.character(diff_pt@data$spi_y15))


saveRDS(diff_pt, "/data/emily/WF/kate/final_data/pt_data_may17.rds") 

#test dataset
test <- readRDS('/data/emily/WF/kate/final_data/pt_data_may17.rds')
coordinates(test)<- c("lon", "lat")
columns<-names(test)
test@proj4string<-diff@crs

spplot(test[,"gw_01"], cex=0.2) 



# data exploration
cv <- unionSpatialPolygons(cv_huc, cv_huc@data$ID)

test <- diff[[13]]
cv <- spTransform(cv, test@crs)

test[test>40] <- NA

test2 <- diff[[13]]
test2 <- crop(test2, cv)


