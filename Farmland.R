# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(maptools)

setwd("C:/ ")

#This script reads in and merges county-extent shapefiles of "farmland of importance" from the 
#California Department of Conservation's Farmland Mapping and Monitring Program (http://www.conservation.ca.gov/dlrp/fmmp/Pages/county_info.aspx).
#Files are posted to their FTP site at ftp://ftp.consrv.ca.gov/pub/dlrp/FMMP/. 
#To use this script download the datafiles of interest, taking note of the file naming convention used below.
#Running the script will produce a state-extent shapefile of California farmland for a single year. 
#Change the year in the file names and repeat for all fmmp years of interest.
#Note that 2014 FMMP GIS data is posted in NAD83, while older GIS data, 1984-2012, remains in NAD 27.

#####################################
### IMPORT COUNTY FMMP SHAPEFILES ###
#####################################

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

###########################
### FIX ATTRIBUTE NAMES ### 
###########################

#Some counties have slightly different names (caplocks) for categories in "data" slot of the spatial dataframe
#This fixes these so data may be merged.

colnames(fmmp_2006.sutter@data) <- c( "polygon_ty" ,"upd_year", "county_nam", "polygon_ac", "Shape_Leng", "Shape_Area")
colnames(fmmp_2006.yuba@data) <- c( "polygon_ty" ,"upd_year", "county_nam", "polygon_ac", "Shape_Leng", "Shape_Area")

#######################################################
### COMBINE FMMP COUNTY POLYGONS INTO SINGLE OBJECT ###
#######################################################

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


##################################################################
### SELECT AGRICULTURAL AREA (farmland and grazing) FROM FMMP  ###
##################################################################

#See the FMMP metadata (ftp://ftp.consrv.ca.gov/pub/dlrp/FMMP/metadata/html/alameda_meta.htm) or 
# http://www.conservation.ca.gov/dlrp/fmmp/mccu/Pages/map_categories.aspx for category descriptions.
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


#########################################
### IMPORT HUC12 WATERSHED SHAPEFILES ### 
########################################

cv_huc <- readOGR(".","cv_huc")
cv_huc$huc12_id <- as.numeric(as.character(cv_huc$HUC12)) #fix the huc12_id which has duplicates

###################################################
### CHECK PROJECTIONS and RECONCILE DIFFERENCES ###
###################################################

proj4string(fmmp.farmselection) #NAD27
proj4string(cv_huc) #NAD83
fmmp.farmselection <- spTransform(fmmp.farmselection,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set to NAD83

##########################################
#WRITE FARMLAND SELECTIONS to SHAPEFILE ###
###########################################
writeOGR(fmmp.farmselection,".","fmmp.farmselection.2006", driver="ESRI Shapefile") #write out a state-extents farmland shapefile


