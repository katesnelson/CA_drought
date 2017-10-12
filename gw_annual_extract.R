library(spacetime)
library(raster)

#This script calculates the difference in predicted groundwater elevations between January of year 1 and JAnueary of year 0. This measure of groundwater was
#elminated from models in favor of a measure of depth to groundwater in January.


#load krigged sto
kd <- readRDS('/data/emily/WF/kate/final_data/krigged_data.rds')

d <- as.data.frame(read.csv('/data/emily/WF/kate/final_data/pt_data_May_05.txt', sep=','), header=T, stringsAsFactors=F)
coordinates (d)<- c("lon", "lat")
columns<-names(d)
d@proj4string<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

#subset kd for first of month
jan <- kd[ ,seq(1,length(kd@time),12)]
jan15 <- kd[,180]

y01 <- raster(jan[,2]) - raster(jan[,1])
y02 <- raster(jan[,3]) - raster(jan[,2])
y03 <- raster(jan[,4]) - raster(jan[,3])
y04 <- raster(jan[,5]) - raster(jan[,4])
y05 <- raster(jan[,6]) - raster(jan[,5]) #2K
y06 <- raster(jan[,7]) - raster(jan[,6]) #-2K
y07 <- raster(jan[,8]) - raster(jan[,7])
y08 <- raster(jan[,9]) - raster(jan[,8])
y09 <- raster(jan[,10]) - raster(jan[,9]) #2K
y10 <- raster(jan[,11]) - raster(jan[,10])
y11 <- raster(jan[,12]) - raster(jan[,11]) #1K
y12 <- raster(jan[,13]) - raster(jan[,12])
y13 <- raster(jan[,14]) - raster(jan[,13])
y14 <- raster(jan15) - raster(jan[,14])

diff <- brick(y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12,y13,y14)
diff_val <- getValues(diff)
hist(diff_val, 1000, xlim=c(-100,100))

#drop values above 500, below -500
diff[diff > 500] <- NA
diff[diff < -500] <- NA

d_prj <- spTransform(d, diff@crs)
diff <- stack(diff)
diff_pt <- extract(diff, d_prj, method="simple", df=T)
rownames(diff_pt) <- d_prj$POINTID
colnames(diff_pt) <- c("POINTID", "gw_01", "gw_02", "gw_03", "gw_04", "gw_05", "gw_06", "gw_07", "gw_08", "gw_09", "gw_10", "gw_11",
                       "gw_12", "gw_13", "gw_14")
d_prj <- as.data.frame(d_prj, stringsAsFactors=F) 
#sapply(d_prj, class) #see many factors
d_prj$POINTID <- as.numeric(d_prj$POINTID)             

pt_full <- merge(d_prj,diff_pt, by="POINTID", all=T)
saveRDS(pt_full, "/data/emily/WF/kate/final_data/pt_data.rds")

test <- readRDS('/data/emily/WF/kate/final_data/pt_data.rds')

