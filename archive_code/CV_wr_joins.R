# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)

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
#do we want to retain grazing????
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


### IMPORT HUC12 ShAPEFILES ### 
cv_huc <- readOGR(".","cv_huc")
cv_huc$huc12_id <- as.numeric(as.character(cv_huc$HUC12)) #fix the huc12_id which has duplicates

### CHECK PROJECTIONS and RECONCILE DIFFERENCES ###
proj4string(fmmp.farmselection) #NAD27
proj4string(cv_huc) #NAD83

#WRITE FARMLAND SELECTIONS to SHAPEFILE ###
fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set projection to NAD83
#writeOGR(fmmp.farmselection,".","fmmp.farmselection.2006", driver="ESRI Shapefile") #write the farmland shapefile


###Calculate fmmp area per huc ####  
cv_huc <- spTransform(cv_huc,CRS("+proj=longlat +ellps=WGS84") )#set projection to WGS84 for sp
fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +ellps=WGS84") ) #set projection to WGS84 for sp

#fmmp shapefile cleanup
names <- c("fl_typ","Acres","county","year","fmmp_leng","fmmp_area") #rename fmmp attributes to avoid duplicate col names
names(fmmp.farmselection) <- names #associate the new names with the original spatialdataframe
ids <-as.character(row.names(fmmp.farmselection)) #create a column for fmmp ids
fmmp.farmselection$fmmp_id <- ids #add fmmp_id column to fmmp spatial dataframe
keep <-c("fl_typ","Acres", "county", "fmmp_id") #names of columns in cv_huc to keep
fmmp.farmselection <- fmmp.farmselection[, (names(fmmp.farmselection) %in% keep)] #reduce fmmp to critical attributes
fmmp.red <- fmmp.farmselection[,4] #reduce fmmp.farmselection to fmmp_id for joins
fmmp.buffer <- gBuffer(fmmp.red, byid=TRUE, width=0) #get's rid of self-intersection issues

#huc shapefile cleanup
keep <-c("huc12_id","AreaAcres", "AreaSqKm", "diverse", "Name", "HUC12") #names of columns in cv_huc to keep
cv_huc <- cv_huc[, (names(cv_huc) %in% keep)] #reduce cv_huc to critical attributes
huc.red <- cv_huc[, 5] #reduce huc dataset to only huc_id for joins
row.names(huc.red)<-as.character(huc.red$huc12_id)
huc.buffer <- gBuffer(huc.red, byid=TRUE, width=0)

#spatial intersection of fmmp and hucs using raster package
huc.buffer<- spTransform(huc.red,CRS("+proj=longlat +ellps=WGS84") )
fmmp.buffer <- spTransform(fmmp.buffer,CRS("+proj=longlat +ellps=WGS84") )
fmmp.huc.bound<- fmmp.buffer[huc.buffer, ] #take the subset of fmmp polygons that intersect cv hucs
#plot(fmmp.huc.bound)
#fmmp.huc.bound2 <- union( fmmp.huc.bound, huc.buffer)

#first do a union of fmmp with itself to create on large polygon
#then calculate the gArea(gIntersection(huc, fmmp))
fmmp.agg <- gUnaryUnion(fmmp.buffer)
proj(fmmp.agg)<-proj(huc.red)
fmmp.area.huc<-union(huc.buffer, fmmp.agg, byid=TRUE)

#fmmp.huc.bound2$id<-as.character(row.names(fmmp.huc.bound2))
#plot(fmmp.huc.bound2)
#plot(huc.red, add = TRUE, col = "red")
##Union of fmmp areas by huc
fmmp.huc.area <- Area(fmmp.agg.bound2)
fmmp.huc.area<-merge(fmmp.huc.area, fmmp.agg.bound2, by = "row.names")
fmmp.agg.huc<- aggregate(fmmp.huc.area, by = "huc12_id", sums=fmmp.huc.area(fmmp.huc.area(sum,'area')))
#fmmp.agg.huc2 <- gUnaryUnion(fmmp.huc.bound2, id = "huc12_id)


plot(fmmp.agg.huc)

##reintroduce the fmmp and huc attribute data
combine.fmmp.huc <-cbind(fmmp.farmselection,overlay.1) #combine the farmland selection data and the overlay dataframe, 
                                                      #this serves to link the attribute data for each fmmp record with the overlay data
combine.fmmp.huc <-left_join(combine.fmmp.huc, cv_huc@data, by = "huc12_id") #join fmmp data with the huc dataset to get huc attributes for each fmmp
merge.fmmp.huc <-merge(fmmp.red,combine.fmmp.huc, by.y ="fmmp_id", by.x = "row.names") #merge the farmland spatial dataframe with the attribute data 
                                                  #from the above join (by common record name)
farm_HUC <- merge.fmmp.huc[!is.na(merge.fmmp.huc$huc12_id), -1] #remove rows with NAs in the spatialdataframe and the weird extra row name column

##write to a shapefile
#farm_HUC <- spTransform(farm_HUC,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set to WGS84 for mapping
#writeOGR(farm_HUC,".","farm_HUC.2006", driver="ESRI Shapefile")#write to an esri shapefile

### IMPORT WATER RIGHTS DATA  ### 

full.pt.2006<- readOGR(".","full.pt.2006")

###Associate WR Data with HUCs####

