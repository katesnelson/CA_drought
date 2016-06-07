# Required packages
library(dplyr)
#library (gstat)
library(sp)
library(rgdal)
#library(INLA)
library(spdep)
library(taRifx.geo)
library(rgeos)
#library(lattice)
#library(latticeExtra)
library(raster)
library(readxl)

setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/central valley")
#setwd("C:/Users/tuan/Documents/Research")

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
cv_huc$huc12_id <- as.numeric(as.character(cv_huc$HUC12)) #fix the huc12_id which has duplicates
NHD <- readOGR(".","NHDArea") #may not need this (if farmland is intersected by an NHD Area 'Stream River' remove it to ease WR aggregation)

### CHECK PROJECTIONS and RECONCILE DIFFERENCES ###
proj4string(fmmp.farmselection) #NAD27
proj4string(cv_huc) #NAD83

#WRITE FARMLAND SELECTIONS to SHAPEFILE ###
fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set projection to NAD83
#writeOGR(fmmp.farmselection,".","fmmp.farmselection.2006", driver="ESRI Shapefile") #write the farmland shapefile


###ASSOCIATE fmmp POLYGONS WITH huc POLYGONS ####  
cv_huc <- spTransform(cv_huc,CRS("+proj=longlat +ellps=WGS84") )#set projection to WGS84 for rgeos
fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +ellps=WGS84") ) #set projection to WGS84 for rgeos

##fmmp shapefile cleanup
names <- c("fl_typ","Acres","county","year","fmmp_leng","fmmp_area") #rename fmmp attributes to avoid duplicate col names
names(fmmp.farmselection) <- names #associate the new names with the original spatialdataframe
ids <-as.character(row.names(fmmp.farmselection)) #create a column for fmmp ids
fmmp.farmselection$fmmp_id <- ids #add fmmp_id column to fmmp spatial dataframe
keep <-c("fl_typ","Acres", "county", "fmmp_id") #names of columns in cv_huc to keep
fmmp.farmselection <- fmmp.farmselection[, (names(fmmp.farmselection) %in% keep)] #reduce fmmp to critical attributes
fmmp.red <- fmmp.farmselection[,7] #reduce fmmp.farmselection to fmmp_id for joins
buffertest <- gBuffer(fmmp.red, byid=TRUE, width=0) #get's rid of self-intersection issues

##huc shapefile cleanup
keep <-c("huc12_id","AreaAcres", "AreaSqKm", "diverse", "Name", "HUC12") #names of columns in cv_huc to keep
cv_huc <- cv_huc[, (names(cv_huc) %in% keep)] #reduce cv_huc to critical attributes
huc.red <- cv_huc[, 5] #reduce huc dataset to only huc_id for joins

##spatial joins
overlay.1 <- over(buffertest, huc.red, returnList=FALSE, minDimension = 2) #return the attribute data for one cv_huc with areal overlap with each fmmp.farmselection (returns a dataframe)
overlay.1.dup <- overlay.1 #duplicate of overlay.1 for comparison with modified overlay below
overlay.2 <- over(buffertest, huc.red, returnList=TRUE, minDimension = 2) #return the attributes of cv_huc at multiple areal intersection locations with each fmmp.farmselection (returns a list)
buffertest <- spTransform(buffertest,CRS("+proj=longlat +ellps=WGS84") ) #set projection to WGS84 for rgeos
huc.red <- spTransform(huc.red,CRS("+proj=longlat +ellps=WGS84") ) #set projection to WGS84 for rgeos

##calculate areas when a many-to-one join occurs
for (i in 1:nrow(overlay.1)){        #for each record in the short overlay dataframe
  tmp <- length(overlay.2[[i]][[1]])     #calculate the number of hucs that intersect an fmmp record from the long overlay list
  if (tmp>1){                               #if there is more than one huc that intersects a fmmp record
    areas <- numeric(tmp)                   
    for (j in 1:tmp){                       #for each intersection
      areas[j] <- gArea(gIntersection(buffertest[i,],      
                                      huc.red[(huc.red@data$huc12_id == overlay.2[[i]][[1]][j]),]))
    }                                       #calculate the area of the intersection between fmmp records and hucs with
                                            #a huc12_id_2 that matches the huc12_id_2 listed in the long overlay list 
                                            #(make sure 2nd index corresponds to attribute column of interest (huc12_id_2))
    overlay.1[i, "huc12_id"] <- overlay.2[[i]][[1]][which(areas==max(areas))]
  }                                         #replace the huc12_id_2 recorded in the first short overlay with the huc12_id_2 
                                            #of the huc that has the greatest area of overlap with the fmmp record
}

##reintroduce the fmmp and huc attribute data
combine.fmmp.huc <-cbind(fmmp.farmselection,overlay.1) #combine the farmland selection data and the overlay dataframe, 
                                                      #this serves to link the attribute data for each fmmp record with the overlay data
combine.fmmp.huc <-left_join(combine.fmmp.huc, cv_huc@data, by = "huc12_id") #join fmmp data with the huc dataset to get huc attributes for each fmmp
merge.fmmp.huc <-merge(fmmp.red,combine.fmmp.huc, by.y ="fmmp_id", by.x = "row.names") #merge the farmland spatial dataframe with the attribute data 
                                                  #from the above join (by common record name)
farm_HUC <- merge.fmmp.huc[!is.na(merge.fmmp.huc$huc12_id), -1] #remove rows with NAs in the spatialdataframe and the weird extra row name column

##write to a shapefile
#farm_HUC <- spTransform(farm_HUC,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set projection to WGS84 for plotting
#writeOGR(farm_HUC,".","farm_HUC.2006", driver="ESRI Shapefile")#write to an esri shapefile

### IMPORT WATER RIGHTS DATA and join to make more managable ### 

full.pt.2006<- readOGR(".","full.pt.2006")


###POINTs in POLYGON JOIN for WR PODS to FMMP polygons with Huc join info###
farm_HUC <- spTransform(farm_HUC,CRS("+proj=longlat +ellps=WGS84") ) #make sure you're in the same coordibate system for join
full.pt.2006 <- spTransform(full.pt.2006,CRS("+proj=longlat +ellps=WGS84") )
pt.fmmp <- over(full.pt.2006,farm_HUC)
pt.fmmp <-merge(full.pt.2006,pt.fmmp, by = "row.names")
pt.in.fmmp <- pt.fmmp[!is.na(pt.fmmp$huc12_id), ] #keep the subset of wr pts located in fmmp polygons
pt.notin.fmmp <- pt.fmmp[is.na(pt.fmmp$huc12_id), ]
#plot(pt.in.fmmp)
plot(pt.notin.fmmp)

###If Point NOT in FMMP Polygon, Associate to Closest FMMP Polygon if Near### 	
#Proj4js.defs["EPSG:3310"] = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs";
farm_HUC<- spTransform(farm_HUC,CRS("+init=epsg:3310") ) #try setting a projection so can specify buffer in meters
fmmp_buffer<-gBuffer(farm_HUC, byid=TRUE, width=250.0) #buffer meters
#plot(fmmp_buffer)
pt.notin.fmmp<- spTransform(pt.notin.fmmp,CRS("+init=epsg:3310") )
pt.near.fmmp <- over(pt.notin.fmmp,fmmp_buffer)
pt.near.fmmp <-merge(full.pt.2006,pt.near.fmmp, by = "row.names")
pt.near.fmmp <-pt.near.fmmp[!is.na(pt.near.fmmp$huc12_id), ]

###If point not in or near FMMP polygon, Associate to nearby NHD waternetwork junction###



# associate water rights pods with the piece of fmmp that contains them or is nearest to them (within bounds of x meters)
#point.in.polygon()
#buffer around pts?
#unionSpatialPolygons() for dissolve