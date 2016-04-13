# Required packages
library(dplyr)
library (gstat)
library(sp)
library(rgdal)
library(INLA)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(lattice)
library(latticeExtra)
library(raster)

setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/central valley")


### IMPORT fmmp (CA farmland of importance) SHAPEFILES ###

fmmp_2006.yuba <- readOGR("2006_FMMP_shape_files","yuba2006")
fmmp_2006.yolo <- readOGR("2006_FMMP_shape_files","yolo2006")
fmmp_2006.ventura <- readOGR("2006_FMMP_shape_files","ventura2006")
fmmp_2006.tulare <- readOGR("2006_FMMP_shape_files","tulare2006")
fmmp_2006.tehama <- readOGR("2006_FMMP_shape_files","tehama2006")
fmmp_2006.sutter <- readOGR("2006_FMMP_shape_files","sutter2006")
fmmp_2006.stanislaus <- readOGR("2006_FMMP_shape_files","stanislaus2006")
fmmp_2006.sonoma <- readOGR("2006_FMMP_shape_files","sonoma2006")
fmmp_2006.solano <- readOGR("2006_FMMP_shape_files","solano2006")
fmmp_2006.siskiyou <- readOGR("2006_FMMP_shape_files","siskiyou2006")
fmmp_2006.sierravalleyarea <- readOGR("2006_FMMP_shape_files","sierravalleyarea2006")
fmmp_2006.shasta <- readOGR("2006_FMMP_shape_files","shasta2006")
fmmp_2006.santacruz <- readOGR("2006_FMMP_shape_files","santacruz2006")
fmmp_2006.santaclara <- readOGR("2006_FMMP_shape_files","santaclara2006")
fmmp_2006.santabarbara <- readOGR("2006_FMMP_shape_files","santabarbara2006")
fmmp_2006.sanmateo <- readOGR("2006_FMMP_shape_files","sanmateo2006")
fmmp_2006.sanluisobispo <- readOGR("2006_FMMP_shape_files","sanluisobispo2006")
fmmp_2006.sanjoaquin <- readOGR("2006_FMMP_shape_files","sanjoaquin2006")
fmmp_2006.sandiego <- readOGR("2006_FMMP_shape_files","sandiego2006")
fmmp_2006.sanbernardino <- readOGR("2006_FMMP_shape_files","sanbernardino2006")
fmmp_2006.sanbenito <- readOGR("2006_FMMP_shape_files","sanbenito2006")
fmmp_2006.sacramento <- readOGR("2006_FMMP_shape_files","sacramento2006")
fmmp_2006.riverside <- readOGR("2006_FMMP_shape_files","riverside2006")
fmmp_2006.placer <- readOGR("2006_FMMP_shape_files","placer2006")
fmmp_2006.orange <- readOGR("2006_FMMP_shape_files","orange2006")
fmmp_2006.nevada <- readOGR("2006_FMMP_shape_files","nevada2006")
fmmp_2006.napa <- readOGR("2006_FMMP_shape_files","napa2006")
fmmp_2006.monterey <- readOGR("2006_FMMP_shape_files","monterey2006")
fmmp_2006.modoc <- readOGR("2006_FMMP_shape_files","modoc2006")
fmmp_2006.merced <- readOGR("2006_FMMP_shape_files","merced2006")
fmmp_2006.mendocino <- readOGR("2006_FMMP_shape_files","mendocino2006")
fmmp_2006.mariposa <- readOGR("2006_FMMP_shape_files","mariposa2006")
fmmp_2006.marin <- readOGR("2006_FMMP_shape_files","marin2006")
fmmp_2006.madera <- readOGR("2006_FMMP_shape_files","madera2006")
fmmp_2006.losangeles <- readOGR("2006_FMMP_shape_files","losangeles2006")
fmmp_2006.lake <- readOGR("2006_FMMP_shape_files","lake2006")
fmmp_2006.kings <- readOGR("2006_FMMP_shape_files","kings2006")
fmmp_2006.kern <- readOGR("2006_FMMP_shape_files","kern2006")
fmmp_2006.imperial <- readOGR("2006_FMMP_shape_files","imperial2006")
fmmp_2006.glenn <- readOGR("2006_FMMP_shape_files","glenn2006")
fmmp_2006.fresno <- readOGR("2006_FMMP_shape_files","fresno2006")
fmmp_2006.eldorado <- readOGR("2006_FMMP_shape_files","eldorado2006")
fmmp_2006.contracosta <- readOGR("2006_FMMP_shape_files","contracosta2006")
fmmp_2006.colusa <- readOGR("2006_FMMP_shape_files","colusa2006")
fmmp_2006.butte <- readOGR("2006_FMMP_shape_files","butte2006")
fmmp_2006.amador <- readOGR("2006_FMMP_shape_files","amador2006")
fmmp_2006.alameda <- readOGR("2006_FMMP_shape_files","alameda2006")

### FIX ATTRIBUTE NAMES ### (slightly different names (caplocks) for categories in slot "data", fix these)

colnames(fmmp_2006.sutter@data) <- c( "polygon_ty" ,"upd_year", "county_nam", "polygon_ac", "Shape_Leng", "Shape_Area")
colnames(fmmp_2006.yuba@data) <- c( "polygon_ty" ,"upd_year", "county_nam", "polygon_ac", "Shape_Leng", "Shape_Area")

### COMBINE fmmp POLYGONS INTO SINGLE OBJECT ###

fmmp_2006.full=rbind.SpatialPolygonsDataFrame(fmmp_2006.alameda,fmmp_2006.amador,fmmp_2006.butte,
                                              fmmp_2006.colusa, fmmp_2006.contracosta, fmmp_2006.eldorado, 
                                              fmmp_2006.fresno, fmmp_2006.glenn, fmmp_2006.imperial,
                                              fmmp_2006.kern, fmmp_2006.kings, fmmp_2006.lake,
                                              fmmp_2006.napa, fmmp_2006.losangeles, fmmp_2006.madera,
                                              fmmp_2006.marin, fmmp_2006.mariposa, fmmp_2006.mendocino,
                                              fmmp_2006.merced, fmmp_2006.modoc, fmmp_2006.monterey,
                                              fmmp_2006.nevada, fmmp_2006.orange, fmmp_2006.placer,
                                              fmmp_2006.riverside, fmmp_2006.sacramento, fmmp_2006.sanbenito,
                                              fmmp_2006.sanbernardino, fmmp_2006.sandiego, fmmp_2006.sanjoaquin,
                                              fmmp_2006.sanluisobispo, fmmp_2006.sanmateo, fmmp_2006.santabarbara,
                                              fmmp_2006.santaclara, fmmp_2006.santacruz, fmmp_2006.shasta,
                                              fmmp_2006.sierravalleyarea, fmmp_2006.siskiyou, fmmp_2006.solano,
                                              fmmp_2006.sonoma, fmmp_2006.stanislaus, fmmp_2006.sutter, 
                                              fmmp_2006.tehama, fmmp_2006.tulare, fmmp_2006.yolo, fmmp_2006.ventura,
                                              fmmp_2006.yuba, fix.duplicated.IDs = TRUE)



### SELECT AGRICULTURAL AREA FROM fmmp (farmland and grazing) ###

fmmp.targets <- c("P", "U", "S", "G", "L", "LP", "I", "N")
fmmp.farmselection <- fmmp_2006.full[fmmp_2006.full@data$polygon_ty%in%fmmp.targets,]

rm(fmmp_2006.alameda,fmmp_2006.amador,fmmp_2006.butte,
   fmmp_2006.colusa, fmmp_2006.contracosta, fmmp_2006.eldorado, 
   fmmp_2006.fresno, fmmp_2006.glenn, fmmp_2006.imperial,
   fmmp_2006.kern, fmmp_2006.kings, fmmp_2006.lake,
   fmmp_2006.napa, fmmp_2006.losangeles, fmmp_2006.madera,
   fmmp_2006.marin, fmmp_2006.mariposa, fmmp_2006.mendocino,
   fmmp_2006.merced, fmmp_2006.modoc, fmmp_2006.monterey,
   fmmp_2006.nevada, fmmp_2006.orange, fmmp_2006.placer,
   fmmp_2006.riverside, fmmp_2006.sacramento, fmmp_2006.sanbenito,
   fmmp_2006.sanbernardino, fmmp_2006.sandiego, fmmp_2006.sanjoaquin,
   fmmp_2006.sanluisobispo, fmmp_2006.sanmateo, fmmp_2006.santabarbara,
   fmmp_2006.santaclara, fmmp_2006.santacruz, fmmp_2006.shasta,
   fmmp_2006.sierravalleyarea, fmmp_2006.siskiyou, fmmp_2006.solano,
   fmmp_2006.sonoma, fmmp_2006.stanislaus, fmmp_2006.sutter, 
   fmmp_2006.tehama, fmmp_2006.tulare, fmmp_2006.yolo, fmmp_2006.ventura,
   fmmp_2006.yuba)


### IMPORT RELEVANT HYDROLOGY ShAPEFILES ### 

cv_huc <- readOGR(".","cv_huc")
# NHD <- readOGR(".","NHDArea") #may not need this (if farmland is intersected by an NHD Area 'Stream River' remove it to ease WR aggregation)

### CHECK PROJECTIONS and RECONCILE DIFFERENCES ###
proj4string(fmmp.farmselection) #NAD27
proj4string(cv_huc) #NAD83

fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set projection to NAD83
writeOGR(fmmp.farmselection,".","fmmp.farmselection.2006", driver="ESRI Shapefile") #write the farmland shapefile



###UPDATE PLOTS#### bleh

#spplot(fmmp.farmselection, col ='black', fill="transparent") +
#layer(cv_huc, c("huc12_id"), col = 'red', fill="transparent")
#spplot(fmmp.farmselection, c("polygon_ac"), col.regions="transparent")
#spplot(cv_huc, c("huc12_id"), col = 'red', col.regions="transparent", add=TRUE)

###ASSOCIATE fmmp POLYGONS WITH huc POLYGONS ####  

overlay.fmmp.huc <- over(fmmp.farmselection, cv_huc, returnList=FALSE, minDimension = 2)#return the attributes of the cv_huc at locations of fmmp.farmselection, return only one huc for each fmmp and base joins on area of intersection
ids <-as.character(rownames(overlay.fmmp.huc)) #create a
names <- c("farmland_type","polygon_ac","county","year","fmmp_leng","fmmp_area") #rename fmmp attributes to avoid duplicate col names
names(fmmp.farmselection) <- names #associate the new names with the original spatialdataframe

buffertest <- gBuffer(fmmp.farmselection, byid=TRUE, width=0)#try a negative buffer of 10-20 meters
overlay.test1 <- over(buffertest, cv_huc, returnList=FALSE, minDimension = 2) #not sure that this did anything

combine.fmmp.huc <-cbind(fmmp.farmselection,overlay.test) #combine the farmland selection data and the overlay dataframe, this serves to link the name of each fmmp record with the overlay data
merge.fmmp.huc <-merge(fmmp.farmselection,combine.fmmp.huc) #merge the farmland selection spatial data with the attribute data from the above merge (by common record name)
farm_HUC <- merge.fmmp.huc[!is.na(merge.fmmp.huc$OBJECTID), ] #remove rows with NAs int he spatialdataframe
writeOGR(farm_HUC,".","farm_HUC_test.2006", driver="ESRI Shapefile")#write to an esri shapefile
plot(farm_HUC)


combine.fmmp.huc <-cbind(fmmp.farmselection,overlay.fmmp.huc) #combine the farmland selection data and the overlay dataframe, this serves to link the name of each fmmp record with the overlay data
merge.fmmp.huc <-merge(fmmp.farmselection,combine.fmmp.huc) #merge the farmland selection spatial data with the attribute data from the above merge (by common record name)
farm_HUC <- merge.fmmp.huc[!is.na(merge.fmmp.huc$OBJECTID), ] #remove rows with NAs int he spatialdataframe
writeOGR(farm_HUC,".","farm_HUC.2006", driver="ESRI Shapefile")#write to an esri shapefile
plot(farm_HUC)

overlay.test2 <- over(buffertest, cv_huc, returnList=TRUE, areaWeighted=TRUE)

overlay.test3 <- over(buffertest, cv_huc, returnList=FALSE, minDimension = 2, fn = max(gArea(gIntersection(buffertest,cv_huc))))

Areaint <- gArea(gIntersection(buffertest,cv_huc, byid=TRUE))
overlay.test3 <- over(buffertest, cv_huc, returnList=FALSE, minDimension = 2, fn = max(Areaint[i]))


cv_huc <- spTransform(cv_huc,CRS("+proj=longlat +ellps=WGS84") )
buffertest <- spTransform(buffertest,CRS("+proj=longlat +ellps=WGS84") ) #set projection to WGS84
overlay.test4 <- spTransform(overlay.test4,CRS("+proj=longlat +ellps=WGS84") ) 
fmmp_to_huc <-data.frame(buffertest@data, HUCid = overlay.test1)
overlay.test4 <- over(buffertest, cv_huc, returnList=TRUE)
for (i in 1:nrow(fmmp_to_huc@polygons)){
  tmp <- length(overlay.test4[[i]])
  if (tmp>1){
    areas <- numeric(tmp)
    for (j in 1:tmp){
      areas[j] <- gArea(gIntersection(buffertest[i,], cv_huc[overlay.test4[[i]][j],]))
    }
    
    #       Next line adds some tiny random jittering because
    #       there is a case (Kawerau) which is an exact tie
    #       in intersection area with two TAs - Rotorua and Whakatane
    
    #areas <- areas * rnorm(tmp,1,0.0001)
    
    fmmp_to_huc[i, "HUCid"] <- overlay.test4[[i]][which(areas==max(areas))]
  }
  
}


# associate water rights pods with the piece of fmmp that contains them or is nearest to them (within bounds of x meters)
#point.in.polygon()
#buffer around pts?
#unionSpatialPolygons() for dissolve