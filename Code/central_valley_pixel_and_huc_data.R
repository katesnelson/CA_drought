# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(maptools)
library(lattice)

setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/central valley")
#setwd("C:/Users/tuan/Documents/Research")

########################
#read in pixel data
##########################
d <- as.data.frame(read.csv('pt_data_May_05.txt', sep=','), header=T, stringsAsFactors=F)  #read in full pixel data set
coordinates (d)<- c("lon", "lat")
columns<-names(d)

#read in central valley huc shapefile
cv_huc <- readOGR(".","cv_huc")
d@proj4string<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") 

#trim pixels to central valley
d_CV <- d[cv_huc, ] #Central valley subset of pixel data
plot(d_CV)

#hist(d_CV$fl_count)
#hist(d_CV$px_tvp_y15)
#hist(d_CV$px_tvp_y6)
#Check classes of columns (spi is factor)

####################################
#read in water right huc level data
######################################
wr <- readOGR(".","full.wrbyhuc2")  #read in full huc dataset
proj4string(wr)
wr<-spTransform(wr,"+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0" )

plot(wr)

wr_dat <-wr@data 
View(wr_dat)#look at the dataframe


#names from creation of the shapefile
names(wr@data) = c("huc12_id", "AreaAcres",  "AreaSqKm",   "HUC12_id",      "Name",   "diverse.06",    "wr_cnt.06",     "rip_cnt.06",    "pre_cnt.06",   
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

huc.wr<-wr #make a copy

#names that we want to keep in final dataset
huc_names.final<- c("huc12_id", "AreaAcres",  "AreaSqKm",    "Name",   "diverse.06",    "wr_cnt.06",     "rip_cnt.06",    "pre_cnt.06",   
                    "ag_cnt.06",     "dom_cnt.06",    "ind_cnt.06",    "fish_cnt.06",   "rec_cnt.06",    "adjud_cnt.06",  "unauth_cnt.06",
                    "diverse.07",    "wr_cnt.07",     "rip_cnt.07",    "pre_cnt.07",   
                    "ag_cnt.07",     "dom_cnt.07",    "ind_cnt.07",    "fish_cnt.07",   "rec_cnt.07",    "adjud_cnt.07",  "unauth_cnt.07",
                    "diverse.08",    "wr_cnt.08",     "rip_cnt.08",    "pre_cnt.08",   
                    "ag_cnt.08",     "dom_cnt.08",    "ind_cnt.08",    "fish_cnt.08",   "rec_cnt.08",    "adjud_cnt.08",  "unauth_cnt.08",
                    "diverse.09",    "wr_cnt.09",     "rip_cnt.09",    "pre_cnt.09",   
                    "ag_cnt.09",     "dom_cnt.09",    "ind_cnt.09",    "fish_cnt.09",   "rec_cnt.09",    "adjud_cnt.09",  "unauth_cnt.09",
                    "diverse.10",    "wr_cnt.10",     "rip_cnt.10",    "pre_cnt.10",   
                    "ag_cnt.10",     "dom_cnt.10",    "ind_cnt.10",    "fish_cnt.10",   "rec_cnt.10",    "adjud_cnt.10",  "unauth_cnt.10",
                    "diverse.11",    "wr_cnt.11",     "rip_cnt.11",    "pre_cnt.11",   
                    "ag_cnt.11",     "dom_cnt.11",    "ind_cnt.11",    "fish_cnt.11",   "rec_cnt.11",    "adjud_cnt.11",  "unauth_cnt.11",
                    "diverse.12",    "wr_cnt.12",     "rip_cnt.12",    "pre_cnt.12",   
                    "ag_cnt.12",     "dom_cnt.12",    "ind_cnt.12",    "fish_cnt.12",   "rec_cnt.12",    "adjud_cnt.12",  "unauth_cnt.12",
                    "diverse.13",    "wr_cnt.13",     "rip_cnt.13",    "pre_cnt.13",   
                    "ag_cnt.13",     "dom_cnt.13",    "ind_cnt.13",    "fish_cnt.13",   "rec_cnt.13",    "adjud_cnt.13",  "unauth_cnt.13",
                    "diverse.14",    "wr_cnt.14",     "rip_cnt.14",    "pre_cnt.14",   
                    "ag_cnt.14",     "dom_cnt.14",    "ind_cnt.14",    "fish_cnt.14",   "rec_cnt.14",    "adjud_cnt.14",  "unauth_cnt.14",
                    "diverse.15",    "wr_cnt.15",     "rip_cnt.15",    "pre_cnt.15",   
                    "ag_cnt.15",     "dom_cnt.15",    "ind_cnt.15",    "fish_cnt.15",   "rec_cnt.15",    "adjud_cnt.15",  "unauth_cnt.15")

#select only data we want to keep, remove duplicate areas, etc...
huc.wr@data <- huc.wr@data[ , (names(huc.wr@data) %in% huc_names.final)]

k<-huc.wr@data 
View(k) #look at the data to confirm selection of correct data and check data classes ---> looks good


#########################
#PLOTS
############################
#plot point and polygon together to confirm correct spatial alignment
plot(huc.wr)
points(d_CV)

#check a few attribute spatial distributions
spplot(huc.wr[ ,"diverse.07"])
spplot(huc.wr[ ,"wr_cnt.06"])
spplot(d_CV [ ,"fl_cd7"])

###############################
#Spatial join of pixel and huc
#################################

