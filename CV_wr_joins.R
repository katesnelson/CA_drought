# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(maptools)

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
#fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set to NAD83
#writeOGR(fmmp.farmselection,".","fmmp.farmselection.2006", driver="ESRI Shapefile") #write the farmland shapefile

#fmmp shapefile cleanup
#fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +ellps=WGS84"))
#names <- c("fl_typ","Acres","county","year","fmmp_leng","fmmp_area") #rename fmmp attributes to avoid duplicate col names
#names(fmmp.farmselection) <- names #associate the new names with the original spatialdataframe
#ids <-as.character(row.names(fmmp.farmselection)) #create a column for fmmp ids
#fmmp.farmselection$fmmp_id <- ids #add fmmp_id column to fmmp spatial dataframe
#keep <-c("fl_typ","Acres", "county", "fmmp_id") #names of columns in cv_huc to keep
#fmmp.farmselection <- fmmp.farmselection[, (names(fmmp.farmselection) %in% keep)] #reduce fmmp to critical attributes
#fmmp.farmselection<- spTransform(fmmp.farmselection,CRS("+init=epsg:3310") ) #try setting a projection so can specify buffer in meters
#fmmp.buffer <- gBuffer(fmmp.farmselection, width=0) #get's rid of self-intersection issues

#huc shapefile cleanup
cv_huc <- spTransform(cv_huc,CRS("+proj=longlat +ellps=WGS84"))
keep <-c("huc12_id","AreaAcres", "AreaSqKm", "diverse", "Name", "HUC12") #names of columns in cv_huc to keep
cv_huc <- cv_huc[, (names(cv_huc) %in% keep)] #reduce cv_huc to critical attributes
#row.names(cv_huc)<-as.character(cv_huc$huc12_id)
#cv_huc<- spTransform(cv_huc, CRS("+init=epsg:3310") ) #try setting a projection so can specify buffer in meters
#huc.buffer <- gBuffer(cv_huc,  width=0, byid=TRUE)

# t<-checkPolygonsHoles(huc.buffer@polygons[[1]])
#t2<-checkPolygonsHoles(fmmp.buffer@polygons[[1]])

#fmmp.area.huc<-raster::union(huc.buffer, fmmp.buffer) 
#for (i in 1:nrow(huc.buffer@data)){        #for each record of huc
#    areas <- -999                  
 #   areas[i] <- gArea(gIntersection(fmmp.buffer,      
 #                                     huc.buffer[i,]))
                                         #calculate the area of the intersection between fmmp and hucs by hucs
#    huc.buffer$area[i] <- areas[i]
#  }                                         #replace the huc12_id_2 recorded in the first short overlay with the huc12_id_2 
  #of the huc that has the greatest area of overlap with the fmmp record


#fmmp.area.huc<-gIntersection(huc.buffer, fmmp.buffer, byid = TRUE)
#test<-fmmp.area.huc@data

#fmmp.area.huc3<-raster::intersect(huc.buffer, fmmp.buffer, dissolve=FALSE)

#plot(fmmp.area.huc)

### IMPORT WATER RIGHTS DATA  ### 

full.pt.2015<- readOGR(".","full.pt.2015")
full.pt.2015 <- spTransform(full.pt.2015,CRS("+proj=longlat +ellps=WGS84"))

###Associate WR Data with HUCs####

pt.huc.2015 <- raster::intersect(full.pt.2015, cv_huc ) #spatial join of POD pts and huc shapefile
pt.huc.dat <- pt.huc.2015@data #data only from join
cv_huc_dup <-cv_huc #duplicate so we don't screw up the original


###for each huc calculate aggregate values of interest ###

#set of placeholders
wr_cnt<-0
rip_cnt <-0
pre_cnt<-0
adjud_cnt<-0
unauth_cnt<-0
ag_cnt<-0
dom_cnt<-0
ind_cnt<-0
fish_cnt<-0
rec_cnt<-0
for (i in 1:nrow(cv_huc_dup)){
  subset<- dplyr::filter(pt.huc.dat, huc12_id == cv_huc_dup$huc12_id[i])
  wr_cnt[i]<- length(subset$huc12_id)
  cv_huc_dup$wr_cnt[i] <- wr_cnt[i]  #add count of water rights for each huc
  rip <- dplyr::filter(subset, !is.na(Riparin))
  rip_cnt[i]<- length(rip$Riparin)
  cv_huc_dup$rip_cnt[i] <- rip_cnt[i] #add count of riparian wrs for each huc
  pre <- dplyr::filter(subset, !is.na(Pre1914))
  pre_cnt[i]<- length(pre$Pre1914)
  cv_huc_dup$pre_cnt[i] <- pre_cnt[i] #add count of pre1914 wrs for each huc
  ag <- dplyr::filter(subset, use == "Agriculture" )
  ag_cnt[i]<- length(ag$use)
  cv_huc_dup$ag_cnt[i] <- ag_cnt[i] #add count of agricultural wrs for each huc
  dom <- dplyr::filter(subset, use == "Domestic" )
  dom_cnt[i]<- length(dom$use)
  cv_huc_dup$dom_cnt[i] <- dom_cnt[i] #add count of domestic wrs for each huc
  ind <- dplyr::filter(subset, use == "Industrial" )
  ind_cnt[i]<- length(ind$use)
  cv_huc_dup$ind_cnt[i] <- ind_cnt[i] #add count of industrial wrs for each huc
  fish <- dplyr::filter(subset, use == "Fish & Wildlife" )
  fish_cnt[i]<- length(fish$use)
  cv_huc_dup$fish_cnt[i] <- fish_cnt[i] #add count of fish & wildlife wrs for each huc
  rec <- dplyr::filter(subset, use == "Recreation & Other" )
  rec_cnt[i]<- length(rec$use)
  cv_huc_dup$rec_cnt[i] <- rec_cnt[i] #add count of recreational wrs for each huc
  adjud <- dplyr::filter(subset, !is.na(CortDcr))
  adjud_cnt[i]<- length(adjud$CortDcr)
  cv_huc_dup$adjud_cnt[i] <- adjud_cnt[i] #add count of court ajudicated wrs for each huc
  unauth <- dplyr::filter(subset, !is.na(Un.Athr))
  unauth_cnt[i]<- length(unauth$Un.Athr)
  cv_huc_dup$unauth_cnt[i] <- unauth_cnt[i] #add count of court ajudicated wrs for each huc
}


cv_huc_dup <- spTransform(cv_huc_dup,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set to NAD83
writeOGR(cv_huc_dup,".","wrbyhuc.2015", driver="ESRI Shapefile") #write the farmland shapefile

###READ in and MERGE all of the Water Right POD Aggregated to HUC Info ###
wrbyhuc.2006<- readOGR(".","wrbyhuc.2006")
wrbyhuc.2007<- readOGR(".","wrbyhuc.2007")
wrbyhuc.2008<- readOGR(".","wrbyhuc.2008")
wrbyhuc.2009<- readOGR(".","wrbyhuc.2009")
wrbyhuc.2010<- readOGR(".","wrbyhuc.2010")
wrbyhuc.2011<- readOGR(".","wrbyhuc.2011")
wrbyhuc.2012<- readOGR(".","wrbyhuc.2012")
wrbyhuc.2013<- readOGR(".","wrbyhuc.2013")
wrbyhuc.2014<- readOGR(".","wrbyhuc.2014")
wrbyhuc.2015<- readOGR(".","wrbyhuc.2015")

full.wrbyhuc <- sp::merge(wrbyhuc.2006, wrbyhuc.2007@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2008@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2009@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2010@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2011@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2012@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2013@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2014@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2015@data, by = "huc12_id")

names(full.wrbyhuc) = c("AreaAcres.06",  "AreaSqKm.06",   "HUC12_id",      "Name.06",       "huc12_id.06",   "diverse.06",    "wr_cnt.06",     "rip_cnt.06",    "pre_cnt.06",   
                       "ag_cnt.06",     "dom_cnt.06",    "ind_cnt.06",    "fish_cnt.06",   "rec_cnt.06",    "adjud_cnt.06",  "unauth_cnt.06",
                       "AreaAcres.07",  "AreaSqKm.07",    "HUC12_id.07",    "Name.07",       "diverse.07",    "wr_cnt.07",     "rip_cnt.07",    "pre_cnt.07",   
                       "ag_cnt.07",     "dom_cnt.07",    "ind_cnt.07",    "fish_cnt.07",   "rec_cnt.07",    "adjud_cnt.07",  "unauth_cnt.07",
                       "AreaAcres.08",  "AreaSqKm.08",     "HUC12_id.08",    "Name.08",         "diverse.08",    "wr_cnt.08",     "rip_cnt.08",    "pre_cnt.08",   
                       "ag_cnt.08",     "dom_cnt.08",    "ind_cnt.08",    "fish_cnt.08",   "rec_cnt.08",    "adjud_cnt.08",  "unauth_cnt.08",
                       "AreaAcres.09",  "AreaSqKm.09",     "HUC12_id.09", "Name.09",        "diverse.09",    "wr_cnt.09",     "rip_cnt.09",    "pre_cnt.09",   
                       "ag_cnt.09",     "dom_cnt.09",    "ind_cnt.09",    "fish_cnt.09",   "rec_cnt.09",    "adjud_cnt.09",  "unauth_cnt.09",
                       "AreaAcres.10",  "AreaSqKm.10",      "HUC12_id.10",  "Name.10",         "diverse.10",    "wr_cnt.10",     "rip_cnt.10",    "pre_cnt.10",   
                       "ag_cnt.10",     "dom_cnt.10",    "ind_cnt.10",    "fish_cnt.10",   "rec_cnt.10",    "adjud_cnt.10",  "unauth_cnt.10",
                       "AreaAcres.11",  "AreaSqKm.11",      "HUC12_id.11",   "Name.11",        "diverse.11",    "wr_cnt.11",     "rip_cnt.11",    "pre_cnt.11",   
                       "ag_cnt.11",     "dom_cnt.11",    "ind_cnt.11",    "fish_cnt.11",   "rec_cnt.11",    "adjud_cnt.11",  "unauth_cnt.11",
                       "AreaAcres.12",  "AreaSqKm.12",      "HUC12_id.12",  "Name.12",        "diverse.12",    "wr_cnt.12",     "rip_cnt.12",    "pre_cnt.12",   
                       "ag_cnt.12",     "dom_cnt.12",    "ind_cnt.12",    "fish_cnt.12",   "rec_cnt.12",    "adjud_cnt.12",  "unauth_cnt.12",
                       "AreaAcres.13",  "AreaSqKm.13",      "HUC12_id.13", "Name.13",         "diverse.13",    "wr_cnt.13",     "rip_cnt.13",    "pre_cnt.13",   
                       "ag_cnt.13",     "dom_cnt.13",    "ind_cnt.13",    "fish_cnt.13",   "rec_cnt.13",    "adjud_cnt.13",  "unauth_cnt.13",
                       "AreaAcres.14",  "AreaSqKm.14",      "HUC12_id.14", "Name.14",          "diverse.14",    "wr_cnt.14",     "rip_cnt.14",    "pre_cnt.14",   
                       "ag_cnt.14",     "dom_cnt.14",    "ind_cnt.14",    "fish_cnt.14",   "rec_cnt.14",    "adjud_cnt.14",  "unauth_cnt.14",
                       "AreaAcres.15",  "AreaSqKm.15",      "HUC12_id.15",  "Name.15",         "diverse.15",    "wr_cnt.15",     "rip_cnt.15",    "pre_cnt.15",   
                       "ag_cnt.15",     "dom_cnt.15",    "ind_cnt.15",    "fish_cnt.15",   "rec_cnt.15",    "adjud_cnt.15",  "unauth_cnt.15")

writeOGR(full.wrbyhuc,".","full.wrbyhuc", driver="ESRI Shapefile") #write the farmland shapefile

###Associate Emily's Pts file with HUCS to get HUCID at each pixel ####

###Associate fmmp shapefile with pixels ###
#Use this to get farmland type and farmland flag for each pixel
# also aggregate(count function) pixels with positive farmland falg by hucs to get approximate area of farmland per huc

#note: we are doing this second step as gArea(gIntersection(of HUCS and fmmps)) to get intersections and associated areas 
#is not working due to a hole in one of the shapefiles. Does not work in ArcGIS either. Can check using the previous script 
#I wrote for identifying hucs in which most of an fmmp shape that spanned several hucs were located.
