library(dplyr)
library(gstat)
library(spacetime) 
library(sp)  
library(xts) 
library(rgdal)
library(raster)
library(Metrics)
library(ggplot2)
library(lattice)
library(rgeos)
library(tidyr)

setwd("C:/")

#This script reads in the pixel data and watershed data for all years and combines them into a single long 
#format dataset to be used in R-INLA models.

########################
###READ IN PIXEL DATA###
########################

df <- readRDS('pt_data_may17.rds')

#CHECK THE DATA
#View(df)
hist(df$fl_cd7)
hist(df$gw_07)
!is.na(any(df[,"gw_14"]))
is.na(all(df[,"gw_14"]))
pix <- df #make a copy
plot(pix)

#Clip pixel dataset to Cental Valley extents 
cv <- shapefile('cv_huc.shp')
#cv <- shapefile('CV_alluvial_bnd.shp')
par("mar")
par(mar=c(1,1,1,1))
plot(cv)
cv <- spTransform(cv, CRS(proj4string(pix))) #set to same proj
proj4string(cv)==proj4string(pix) #double check proj
pix <- pix[cv,] #Clip
plot(pix)
spplot(pix[!is.na(pix$gw_12), "gw_12"], cex=0.5)

################################################
###READ IN WATER RIGHTS, WATERSHED LEVEL DATA###
################################################

#wr <- readOGR('/home/kate/CA/data/',"full.wrbyhuc2")  #read in full huc dataset
wr <- readRDS('finalfullwrbyhuc.rds')  #read in full huc dataset
proj4string(wr)
wr<-spTransform(wr,CRS(proj4string(pix)))

#CHECKS
#plot(wr)
#wr_dat <-wr@data 
#View(wr_dat)#look at the dataframe


#names from creation of the shapefile
names(wr@data) = c("huc12_id", "AreaAcres",  "AreaSqKm",   "HUC12",      "Name",    "diverse.06",    "wr_cnt.06", "approp_cnt.06", "junior_cnt.06",   "rip_cnt.06",    "pre_cnt.06",   
                   "ag_cnt.06", "ag_approp_cnt.06", "ag_rip_cnt.06", "ag_pre_cnt.06", "am_cnt.06", "dom_cnt.06",    "ind_cnt.06",    "fish_cnt.06",   "rec_cnt.06",    "adjud_cnt.06",  "unauth_cnt.06",
                   "AreaAcres.07",  "AreaSqKm.07",    "HUC12_id.07",    "Name.07",       "diverse.07",    "wr_cnt.07",  "approp_cnt.07", "junior_cnt.07",   "rip_cnt.07",    "pre_cnt.07",   
                   "ag_cnt.07", "ag_approp_cnt.07", "ag_rip_cnt.07", "ag_pre_cnt.07",  "am_cnt.07",  "dom_cnt.07",    "ind_cnt.07",    "fish_cnt.07",   "rec_cnt.07",    "adjud_cnt.07",  "unauth_cnt.07",
                   "AreaAcres.08",  "AreaSqKm.08",     "HUC12_id.08",    "Name.08",         "diverse.08",    "wr_cnt.08",  "approp_cnt.08", "junior_cnt.08",   "rip_cnt.08",    "pre_cnt.08",   
                   "ag_cnt.08",  "ag_approp_cnt.08", "ag_rip_cnt.08", "ag_pre_cnt.08",  "am_cnt.08", "dom_cnt.08",    "ind_cnt.08",    "fish_cnt.08",   "rec_cnt.08",    "adjud_cnt.08",  "unauth_cnt.08",
                   "AreaAcres.09",  "AreaSqKm.09",     "HUC12_id.09", "Name.09",        "diverse.09",    "wr_cnt.09", "approp_cnt.09", "junior_cnt.09",    "rip_cnt.09",    "pre_cnt.09",   
                   "ag_cnt.09",  "ag_approp_cnt.09", "ag_rip_cnt.09", "ag_pre_cnt.09", "am_cnt.09",  "dom_cnt.09",    "ind_cnt.09",    "fish_cnt.09",   "rec_cnt.09",    "adjud_cnt.09",  "unauth_cnt.09",
                   "AreaAcres.10",  "AreaSqKm.10",      "HUC12_id.10",  "Name.10",         "diverse.10",    "wr_cnt.10",  "approp_cnt.10", "junior_cnt.10",   "rip_cnt.10",    "pre_cnt.10",   
                   "ag_cnt.10",  "ag_approp_cnt.10", "ag_rip_cnt.10", "ag_pre_cnt.10", "am_cnt.10",  "dom_cnt.10",    "ind_cnt.10",    "fish_cnt.10",   "rec_cnt.10",    "adjud_cnt.10",  "unauth_cnt.10",
                   "AreaAcres.11",  "AreaSqKm.11",      "HUC12_id.11",   "Name.11",        "diverse.11",    "wr_cnt.11",  "approp_cnt.11", "junior_cnt.11",   "rip_cnt.11",    "pre_cnt.11",   
                   "ag_cnt.11",   "ag_approp_cnt.11", "ag_rip_cnt.11", "ag_pre_cnt.11",  "am_cnt.11","dom_cnt.11",    "ind_cnt.11",    "fish_cnt.11",   "rec_cnt.11",    "adjud_cnt.11",  "unauth_cnt.11",
                   "AreaAcres.12",  "AreaSqKm.12",      "HUC12_id.12",  "Name.12",        "diverse.12",    "wr_cnt.12",   "approp_cnt.12", "junior_cnt.12",  "rip_cnt.12",    "pre_cnt.12",   
                   "ag_cnt.12",   "ag_approp_cnt.12", "ag_rip_cnt.12", "ag_pre_cnt.12",  "am_cnt.12","dom_cnt.12",    "ind_cnt.12",    "fish_cnt.12",   "rec_cnt.12",    "adjud_cnt.12",  "unauth_cnt.12",
                   "AreaAcres.13",  "AreaSqKm.13",      "HUC12_id.13", "Name.13",         "diverse.13",    "wr_cnt.13",    "approp_cnt.13", "junior_cnt.13", "rip_cnt.13",    "pre_cnt.13",   
                   "ag_cnt.13", "ag_approp_cnt.13", "ag_rip_cnt.13", "ag_pre_cnt.13",   "am_cnt.13", "dom_cnt.13",    "ind_cnt.13",    "fish_cnt.13",   "rec_cnt.13",    "adjud_cnt.13",  "unauth_cnt.13",
                   "AreaAcres.14",  "AreaSqKm.14",      "HUC12_id.14", "Name.14",          "diverse.14",    "wr_cnt.14",  "approp_cnt.14", "junior_cnt.14",   "rip_cnt.14",    "pre_cnt.14",   
                   "ag_cnt.14",  "ag_approp_cnt.14", "ag_rip_cnt.14", "ag_pre_cnt.14",  "am_cnt.14", "dom_cnt.14",    "ind_cnt.14",    "fish_cnt.14",   "rec_cnt.14",    "adjud_cnt.14",  "unauth_cnt.14",
                   "AreaAcres.15",  "AreaSqKm.15",      "HUC12_id.15",  "Name.15",         "diverse.15",    "wr_cnt.15",   "approp_cnt.15", "junior_cnt.15",  "rip_cnt.15",    "pre_cnt.15",   
                   "ag_cnt.15",  "ag_approp_cnt.15", "ag_rip_cnt.15", "ag_pre_cnt.15",   "am_cnt.15","dom_cnt.15",    "ind_cnt.15",    "fish_cnt.15",   "rec_cnt.15",    "adjud_cnt.15",  "unauth_cnt.15")
huc.wr<-wr #make a copy

#names that we want to keep in final dataset (remove duplicate areas, hucids, etc...)
huc_names.final<- c("huc12_id", "AreaAcres",  "AreaSqKm",   "HUC12",      "Name",    
                    "diverse.06",    "wr_cnt.06", "approp_cnt.06", "junior_cnt.06",   "rip_cnt.06",    "pre_cnt.06",   
                    "ag_cnt.06", "ag_approp_cnt.06", "ag_rip_cnt.06", "ag_pre_cnt.06", "am_cnt.06", "dom_cnt.06",    "ind_cnt.06",    "fish_cnt.06",   "rec_cnt.06",    "adjud_cnt.06",  "unauth_cnt.06",
                    "diverse.07",    "wr_cnt.07",  "approp_cnt.07", "junior_cnt.07",   "rip_cnt.07",    "pre_cnt.07",   
                    "ag_cnt.07", "ag_approp_cnt.07", "ag_rip_cnt.07", "ag_pre_cnt.07",  "am_cnt.07",  "dom_cnt.07",    "ind_cnt.07",    "fish_cnt.07",   "rec_cnt.07",    "adjud_cnt.07",  "unauth_cnt.07",
                    "diverse.08",    "wr_cnt.08",  "approp_cnt.08", "junior_cnt.08",   "rip_cnt.08",    "pre_cnt.08",   
                    "ag_cnt.08",  "ag_approp_cnt.08", "ag_rip_cnt.08", "ag_pre_cnt.08",  "am_cnt.08", "dom_cnt.08",    "ind_cnt.08",    "fish_cnt.08",   "rec_cnt.08",    "adjud_cnt.08",  "unauth_cnt.08",
                    "diverse.09",    "wr_cnt.09", "approp_cnt.09", "junior_cnt.09",    "rip_cnt.09",    "pre_cnt.09",   
                    "ag_cnt.09",  "ag_approp_cnt.09", "ag_rip_cnt.09", "ag_pre_cnt.09", "am_cnt.09",  "dom_cnt.09",    "ind_cnt.09",    "fish_cnt.09",   "rec_cnt.09",    "adjud_cnt.09",  "unauth_cnt.09",
                    "diverse.10",    "wr_cnt.10",  "approp_cnt.10", "junior_cnt.10",   "rip_cnt.10",    "pre_cnt.10",   
                    "ag_cnt.10",  "ag_approp_cnt.10", "ag_rip_cnt.10", "ag_pre_cnt.10", "am_cnt.10",  "dom_cnt.10",    "ind_cnt.10",    "fish_cnt.10",   "rec_cnt.10",    "adjud_cnt.10",  "unauth_cnt.10",
                    "diverse.11",    "wr_cnt.11",  "approp_cnt.11", "junior_cnt.11",   "rip_cnt.11",    "pre_cnt.11",   
                    "ag_cnt.11",   "ag_approp_cnt.11", "ag_rip_cnt.11", "ag_pre_cnt.11",  "am_cnt.11","dom_cnt.11",    "ind_cnt.11",    "fish_cnt.11",   "rec_cnt.11",    "adjud_cnt.11",  "unauth_cnt.11",
                    "diverse.12",    "wr_cnt.12",   "approp_cnt.12", "junior_cnt.12",  "rip_cnt.12",    "pre_cnt.12",   
                    "ag_cnt.12",   "ag_approp_cnt.12", "ag_rip_cnt.12", "ag_pre_cnt.12",  "am_cnt.12","dom_cnt.12",    "ind_cnt.12",    "fish_cnt.12",   "rec_cnt.12",    "adjud_cnt.12",  "unauth_cnt.12",
                    "diverse.13",    "wr_cnt.13",    "approp_cnt.13", "junior_cnt.13", "rip_cnt.13",    "pre_cnt.13",   
                    "ag_cnt.13", "ag_approp_cnt.13", "ag_rip_cnt.13", "ag_pre_cnt.13",   "am_cnt.13", "dom_cnt.13",    "ind_cnt.13",    "fish_cnt.13",   "rec_cnt.13",    "adjud_cnt.13",  "unauth_cnt.13",
                    "diverse.14",    "wr_cnt.14",  "approp_cnt.14", "junior_cnt.14",   "rip_cnt.14",    "pre_cnt.14",   
                    "ag_cnt.14",  "ag_approp_cnt.14", "ag_rip_cnt.14", "ag_pre_cnt.14",  "am_cnt.14", "dom_cnt.14",    "ind_cnt.14",    "fish_cnt.14",   "rec_cnt.14",    "adjud_cnt.14",  "unauth_cnt.14",
                    "diverse.15",    "wr_cnt.15",   "approp_cnt.15", "junior_cnt.15",  "rip_cnt.15",    "pre_cnt.15",   
                    "ag_cnt.15",  "ag_approp_cnt.15", "ag_rip_cnt.15", "ag_pre_cnt.15",   "am_cnt.15","dom_cnt.15",    "ind_cnt.15",    "fish_cnt.15",   "rec_cnt.15",    "adjud_cnt.15",  "unauth_cnt.15")


#select only data we want to keep, remove duplicate areas, etc...
huc.wr@data <- huc.wr@data[ , (names(huc.wr@data) %in% huc_names.final)]

#checks
#k<-huc.wr@data 
#View(k) #look at the data to confirm selection of correct data and check data classes 

#############################
#Fixing the AreaSqKm Values#
############################
huc.wr@data$AreaSqKm<-area(huc.wr) #in square meters
huc.wr@data$AreaSqKm<-huc.wr@data$AreaSqKm/1000000 #convert to sq km

#########################
#PLOTS
############################
#plot point and polygon together to confirm correct spatial alignment
plot(huc.wr)
points(pix, cex=0.05, col="blue")

#check a few attribute spatial distributions for both resolutions

spplot(pix [ ,"fl_cd7"])
spplot(pix [!is.na(pix$gw_11) ,"gw_11"])
spplot(huc.wr[ ,"wr_cnt.06"])


###############################
#Spatial join of pixel and huc
#################################

#full.dat<- over(pix, huc.wr, returnList=F) #sp::over joins the attribute data but returns a dataframe not SPDF, which is annoying...
full.dat<- intersect(pix, huc.wr) #spatial join of huc and pixel data
test <- full.dat #copy
spplot(test[ , "wr_cnt.06"], cex=0.2) #testing to make sure huc data properly joined pixels --> seems ok
spplot(test[ , "ag_lulc7"], cex=0.2) #testing to make sure huc data properly joined pixels --> seems ok
spplot(test[!is.na(test$gw_12), "gw_12"], cex=0.2) # testing the groundwater info


#############################################
#Build long data format
##############################################
View(full.dat@data)

#select subset of desired data for which data is available at all years 2007-2014
check<-full.dat@data[ , c(-4,-(13:16), -25, -(26:30),-32, -33, -42, -52, -61, -70, -(72:77), -(90:107), -(244:260))]

#reorder columns to consist format (ascending year grouped by variable)
sorttest<-check[ ,c(1,2,3,63:66,11,10,9,8,7,6,5,4,12:62,
                    67,84,101,118,135,152,169,186,
                    68,85,102,119,136,153,170,187,
                    69,86,103,120,137,154,171,188,
                    70,87,104,121,138,155,172,189,
                    71,88,105,122,139,156,173,190,
                    72,89,106,123,140,157,174,191,
                    73,90,107,124,141,158,175,192,
                    74,91,108,125,142,159,176,193,
                    75,92,109,126,143,160,177,194,
                    76,93,110,127,144,161,178,195,
                    77,94,111,128,145,162,179,196,
                    78,95,112,129,146,163,180,197,
                    79,96,113,130,147,164,181,198,
                    80,97,114,131,148,165,182,199,
                    81,98,115,132,149,166,183,200,
                    82,99,116,133,150,167,184,201,
                    83,100,117,134,151,168,185,202)]

sorttest<-sorttest[ ,c(1:7,24,33,58,8:23,25:32,34:57,59:202)]

###convert from wide to long format

#if you rerun gather on the longtest1 it multiplies the number of rows, no easy option for mutiple gathers
#reshape is another option, but requires consistent and recognizable naming conventions
#http://stackoverflow.com/questions/25925556/gather-multiple-sets-of-columns-with-tidyr 

#Do several independent gathers then join based on pointid and date, which I'll add before joining

longtest1<-tidyr::gather(sorttest[ ,c(1:10,11:18)], "LUYR", "lu", c(11:18)) #this gathers properly, combine lu cols and retain identifiers
longtest1<-mutate(longtest1, year= as.numeric(substr(longtest1$LUYR, 3,4))) #add a year column based on column names

longtest2<-tidyr::gather(sorttest[ ,c(2,19:26)], "SPIYR", "spi", c(2:9)) #this gathers properly,combine spi cols and keep POINTID as identifier
longtest2<-mutate(longtest2, year= as.numeric(substr(longtest2$SPIYR, 6,7))) #add a year column

longtest3<-tidyr::gather(sorttest[ ,c(2,27:34)], "pxtvpYR", "px_tvp", c(2:9)) #this gathers properly
longtest3<-mutate(longtest3, year= as.numeric(substr(longtest3$pxtvpYR, 9, 10))) #add a year column

longtest4<-tidyr::gather(sorttest[ ,c(2,35:42)], "flcdyr", "fl_cd", c(2:9)) #this gathers properly
longtest4<-mutate(longtest4, year= as.numeric(substr(longtest4$flcdyr, 6,7))) #add a year column

longtest5<-tidyr::gather(sorttest[ ,c(2,43:50)], "cdyr", "cd", c(2:9)) #this gathers properly
longtest5<-mutate(longtest5, year= as.numeric(substr(longtest5$cdyr, 3,4))) #add a year column

longtest6<-tidyr::gather(sorttest[ ,c(2,51:58)], "aglulcyr", "aglulc", c(2:9)) #this gathers properly
longtest6<-mutate(longtest6, year= as.numeric(substr(longtest6$aglulcyr, 8,9))) #add a year column

longtest7<-tidyr::gather(sorttest[ ,c(2,59:66)], "gwyr", "gw", c(2:9)) #this gathers properly
longtest7<-mutate(longtest7, year= as.numeric(substr(longtest7$gwyr, 4,5))) #add a year column

longtest8<-tidyr::gather(sorttest[ ,c(2,67:74)], "diverseyr", "diverse", c(2:9)) #this gathers properly
longtest8<-mutate(longtest8, year= as.numeric(substr(longtest8$diverseyr, 9,10))) #add a year column

longtest9<-tidyr::gather(sorttest[ ,c(2,75:82)], "wrcntyr", "wr_cnt", c(2:9)) #this gathers properly
longtest9<-mutate(longtest9, year= as.numeric(substr(longtest9$wrcntyr, 8,9))) #add a year column

longtest10<-tidyr::gather(sorttest[ ,c(2,83:90)], "appropcntyr", "approp_cnt", c(2:9)) #this gathers properly
longtest10<-mutate(longtest10, year= as.numeric(substr(longtest10$appropcntyr, 12,13))) #add a year column

longtest11<-tidyr::gather(sorttest[ ,c(2,91:98)], "juniorcntyr", "junior_cnt", c(2:9)) #this gathers properly
longtest11<-mutate(longtest11, year= as.numeric(substr(longtest11$juniorcntyr, 12,13))) #add a year column

longtest12<-tidyr::gather(sorttest[ ,c(2,99:106)], "ripcntyr", "rip_cnt", c(2:9)) #this gathers properly
longtest12<-mutate(longtest12, year= as.numeric(substr(longtest12$ripcntyr, 9,10))) #add a year column

longtest13<-tidyr::gather(sorttest[ ,c(2,107:114)], "precntyr", "pre_cnt", c(2:9)) #this gathers properly
longtest13<-mutate(longtest13, year= as.numeric(substr(longtest13$precntyr, 9,10))) #add a year column

longtest14<-tidyr::gather(sorttest[ ,c(2,115:122)], "agcntyr", "ag_cnt", c(2:9)) #this gathers properly
longtest14<-mutate(longtest14, year= as.numeric(substr(longtest14$agcntyr, 8,9))) #add a year column

longtest15<-tidyr::gather(sorttest[ ,c(2,123:130)], "agappropcntyr", "agapprop_cnt", c(2:9)) #this gathers properly
longtest15<-mutate(longtest15, year= as.numeric(substr(longtest15$agappropcntyr, 15,16))) #add a year column

longtest16<-tidyr::gather(sorttest[ ,c(2,131:138)], "agripcntyr", "agrip_cnt", c(2:9)) #this gathers properly
longtest16<-mutate(longtest16, year= as.numeric(substr(longtest16$agripcntyr, 12,13))) #add a year column

longtest17<-tidyr::gather(sorttest[ ,c(2,139:146)], "agprecntyr", "agpre_cnt", c(2:9)) #this gathers properly
longtest17<-mutate(longtest17, year= as.numeric(substr(longtest17$agprecntyr, 12,13))) #add a year column

longtest18<-tidyr::gather(sorttest[ ,c(2,147:154)], "amcntyr", "am_cnt", c(2:9)) #this gathers properly
longtest18<-mutate(longtest18, year= as.numeric(substr(longtest18$amcntyr, 8,9))) #add a year column

longtest19<-tidyr::gather(sorttest[ ,c(2,155:162)], "domcntyr", "dom_cnt", c(2:9)) #this gathers properly
longtest19<-mutate(longtest19, year= as.numeric(substr(longtest19$domcntyr, 9,10))) #add a year column

longtest20<-tidyr::gather(sorttest[ ,c(2,163:170)], "indcntyr", "ind_cnt", c(2:9)) #this gathers properly
longtest20<-mutate(longtest20, year= as.numeric(substr(longtest20$indcntyr, 9,10))) #add a year column

longtest21<-tidyr::gather(sorttest[ ,c(2,171:178)], "fishcntyr", "fish_cnt", c(2:9)) #this gathers properly
longtest21<-mutate(longtest21, year= as.numeric(substr(longtest21$fishcntyr, 10,11))) #add a year column

longtest22<-tidyr::gather(sorttest[ ,c(2,179:186)], "reccntyr", "rec_cnt", c(2:9)) #this gathers properly
longtest22<-mutate(longtest22, year= as.numeric(substr(longtest22$reccntyr, 9,10))) #add a year column

longtest23<-tidyr::gather(sorttest[ ,c(2,187:194)], "adjudcntyr", "adjud_cnt", c(2:9)) #this gathers properly
longtest23<-mutate(longtest23, year= as.numeric(substr(longtest23$adjudcntyr, 11,12))) #add a year column

longtest24<-tidyr::gather(sorttest[ ,c(2,195:202)], "unauthcntyr", "unauth_cnt", c(2:9)) #this gathers properly
longtest24<-mutate(longtest24, year= as.numeric(substr(longtest24$unauthcntyr, 12,13))) #add a year column

#
#join long format variables (could use bind_cols but want to ensure that all location-time pairs match correctly)
full.long<- left_join(longtest1,longtest2, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest3, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest4, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest5, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest6, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest7, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest8, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest9, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest10, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest11, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest12, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest13, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest14, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest15, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest16, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest17, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest18, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest19, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest20, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest21, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest22, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest23, by=c("POINTID","year"))
full.long<- left_join(full.long,longtest24, by=c("POINTID","year"))

#remove unneccessary columns
full.dat.long<-full.long[ ,c(-11,-14, -16,-18,-20,-22,-24,-26,-28,-30,-32,-34,-36,-38,-40,-42,-44,-46,-48,-50,-52,-54,-56,-58)]
#rearrange columns
full.dat.long<-full.dat.long[ ,c(1:10,12,11,13:35)]

#########################################################
#Add densities and percentages and correct data classes##
#########################################################

#pixel res 1 km for calculating farmland area
full.dat.long <- mutate(full.dat.long, fl_areaKm = fl_count*1, wr_dens = wr_cnt/AreaSqKm, ag_dens = ag_cnt/fl_areaKm, 
                        rip_perc = rip_cnt/wr_cnt*100, pre_perc = pre_cnt/wr_cnt*100, approp_perc = approp_cnt/wr_cnt*100,
                        agrip_perc = agrip_cnt/ag_cnt*100, agpre_perc = agpre_cnt/ag_cnt*100, agapprop_perc = agapprop_cnt/ag_cnt*100,
                        ag_perc = ag_cnt/wr_cnt*100, am_perc = am_cnt/wr_cnt*100, dom_perc = dom_cnt/wr_cnt*100, ind_perc = ind_cnt/wr_cnt*100, 
                        fish_perc=fish_cnt/wr_cnt*100, rec_perc=rec_cnt/wr_cnt*100, adjud_perc = adjud_cnt/wr_cnt*100, unauth_perc=unauth_cnt/wr_cnt*100)

full.dat.long$aglulc<-as.factor(full.dat.long$aglulc) #landuse categories in aglulc are: "Barren & Fallow", "Grasses", "Grains", "Row Crops","Fruits and Nuts","Uncultivated Cover" (respectively)

################################
#Full long data as a Shapefile##
################################

cv.dat<-full.dat
cv.dat<- cv.dat[ ,2] #retain locations and POINTID only
cv.dat<-merge(cv.dat,full.dat.long, by = "POINTID", duplicateGeom=TRUE)

#plot to check joins
plot(cv.dat, cex = 0.005)

##########################################
#Write it out
##########################################
#writeOGR(cv.dat, "/home/kate/CA/data", "cv.dat05312017", driver = "ESRI Shapefile")
saveRDS(cv.dat, file = "/home/kate/CA/data/cv.dat06072017.rds")
cv.dat<-readRDS('cv.dat06072017.rds')

#######################################
#Plots checks and Basic Stat checks
#######################################

sub<-cv.dat
sub<-sub[sub@data$year==07,] #pull out one year at a time
spplot(sub[ ,"AreaSqKm"], cex=0.5)
spplot(sub[ ,"wr_cnt"], cex=0.5)
spplot(sub[ ,"wr_dens"], cex=0.5)
spplot(sub[ ,"pre_cnt"], cex=0.5)
spplot(sub[ , "ag_cnt"], cex=0.5)
spplot(sub[ , "adjud_cnt"], cex=0.5)
spplot(sub[!is.na(sub$gw) , "gw"], cex=0.5) 
spplot(sub[!is.na(sub$agrip_perc) , "agrip_perc"], cex=0.5) 
spplot(sub[!is.na(sub$agpre_perc) , "agpre_perc"], cex=0.5) 
spplot(sub[ ,"fl_cd"])
spplot(sub[ ,"adjud_perc"])
spplot(sub[ ,"pre_perc"])
spplot(sub[ ,"px_tvp"])
spplot(sub[, "aglulc"])


#basic stats and quality tests on dataset
any(is.na(sub@data$AreaSqKm))
hist(sub@data$AreaSqKm)
min(sub@data$AreaSqKm) 
max(sub@data$AreaSqKm) 

any(is.na(sub@data$fl_count))
hist(sub@data$fl_count)
min(sub@data$fl_count) #0 --will result in Inf when ag_cnt divided by this 
max(sub@data$fl_count) 

any(is.na(sub@data$spi))
any(is.na(sub@data$px_tvp))
any(is.na(sub@data$fl_cd))
any(is.na(sub@data$aglulc))
any(is.na(sub@data$gw)) #TRUE
any(is.na(sub@data$wr_cnt))
any(is.na(sub@data$rip_cnt))
any(is.na(sub@data$pre_cnt))
any(is.na(sub@data$approp_cnt))
any(is.na(sub@data$ag_cnt))


any(is.na(sub@data$wr_dens))
hist(sub@data$wr_dens)
mean(sub@data$wr_dens)
min(sub@data$wr_dens) #0 --> because some HUCs have no surface water rights
max(sub@data$wr_dens)
any(is.na(sub@data$ag_dens)) #TRUE

hist(sub@data$ag_dens[sub@data$ag_dens!=0])
hist(sub@data$agrip_perc)
hist(sub@data$agpre_perc)
hist(sub@data$agapprop_perc)