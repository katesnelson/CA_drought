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


########################
#read in pixel data
##########################

df <- readRDS('/data/emily/WF/kate/final_data/pt_data.rds')
#fix data classes
#df$spi_y3<- as.numeric(as.character(df$spi_y3)) #spi from factor to numeric
#df$spi_y4<- as.numeric(as.character(df$spi_y4)) #spi from factor to numeric
#df$spi_y5<- as.numeric(as.character(df$spi_y5)) #spi from factor to numeric
#df$spi_y6<- as.numeric(as.character(df$spi_y6)) #spi from factor to numeric
#df$spi_y7<- as.numeric(as.character(df$spi_y7)) #spi from factor to numeric
#df$spi_y8<- as.numeric(as.character(df$spi_y8)) #spi from factor to numeric
#df$spi_y9<- as.numeric(as.character(df$spi_y9)) #spi from factor to numeric
#df$spi_y10<- as.numeric(as.character(df$spi_y10)) #spi from factor to numeric
#df$spi_y11<- as.numeric(as.character(df$spi_y11)) #spi from factor to numeric
#df$spi_y12<- as.numeric(as.character(df$spi_y12)) #spi from factor to numeric
#df$spi_y13<- as.numeric(as.character(df$spi_y13)) #spi from factor to numeric
#df$spi_y14<- as.numeric(as.character(df$spi_y14)) #spi from factor to numeric
#df$spi_y15<- as.numeric(as.character(df$spi_y15)) #spi from factor to numeric

#check the data
#View(df)
hist(df$fl_cd7)
hist(df$gw_07)
!is.na(any(df[,"gw_14"]))
is.na(all(df[,"gw_14"]))
pix <- df #make a copy
plot(pix)

#convert to spatial object --> 082416 now starting with an SPDF instead of text file
#coordinates(pix) <- c("lon", "lat")
#pix@proj4string <- CRS('+init=epsg:3310')
#plot(pix, cex=0.005)
#gridded(pix) = TRUE #error --> suggested tolerance minimum: 0.196956 Error in points2grid(points, tolerance, round) : dimension 1 : coordinate intervals are not constant

#subset point dataset with CV shapefile
cv <- shapefile('/data/emily/WF/kate/shp/cv_huc.shp')
plot(cv)
cv <- spTransform(cv, CRS(proj4string(pix))) #set to same proj
proj4string(cv)==proj4string(pix) #double check proj

pix <- pix[cv,]
plot(pix)

####################################
#read in water right huc level data
######################################
wr <- readOGR('/home/kate/CA/data/',"full.wrbyhuc2")  #read in full huc dataset
proj4string(wr)
wr<-spTransform(wr,CRS(proj4string(pix)))

#checks
#plot(wr)
#wr_dat <-wr@data 
#View(wr_dat)#look at the dataframe


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

#checks
#k<-huc.wr@data 
#View(k) #look at the data to confirm selection of correct data and check data classes ---> looks good

#########################
#PLOTS
############################
#plot point and polygon together to confirm correct spatial alignment
plot(huc.wr)
points(pix, cex=0.005, col="blue")

#check a few attribute spatial distributions for both resolutions
spplot(huc.wr[ ,"diverse.07"])
spplot(pix [ ,"fl_cd7"])
spplot(pix [ ,"gw_11"])
spplot(huc.wr[ ,"wr_cnt.06"])


###############################
#Spatial join of pixel and huc
#################################

#full.dat<- over(pix, huc.wr, returnList=F) #sp::over joins the attribute data but returns a dataframe not SPDF, which is annoying...
full.dat<- intersect(pix, huc.wr) #spatial join of huc and pixel data
test <- full.dat #copy
spplot(test[ , "wr_cnt.06"], cex=0.2) #testing to make sure huc data properly joined pixels --> seems ok
spplot(test[ , "ag_lulc7"], cex=0.2) #testing to make sure huc data properly joined pixels --> seems ok
spplot(test[ , "gw_12"], cex=0.2) # testing the groundwater info --> looks really weird --> check this out!!! --> FIXED 082416

#gridded(test) = TRUE
#g<- points2grid(test, tolerance=0.00820492, round=NULL)

#############################################
#Build long data format
##############################################
View(full.dat@data)

#select subset of data for which all variables are available at all years
check<-full.dat@data[ , c(-4,-(13:16), -25, -(26:30),-32, -33, -42, -52, -61, -70, -(72:77), -(90:100))]
check<-check[ , -(155:165)]

#reorder columns to consist format (ascending year grouped by variable)
sorttest<-check[ ,c(1,2,3,63:66,11,10,9,8,7,6,5,4,12,13:62,
                    67,78,89,100,111,122,133,144,
                    68,79,90,101,112,123,134,145,
                    69,80,91,102,113,124,135,146,
                    70,81,92,103,114,125,136,147,
                    71,82,93,104,115,126,137,148,
                    72,83,94,105,116,127,138,149,
                    73,84,95,106,117,128,139,150,
                    74,85,96,107,118,129,140,151,
                    75,86,97,108,119,130,141,152,
                    76,87,98,109,120,131,142,153,
                    77,88,99,110,121,132,143,154 )]

sorttest<-sorttest[ ,c(1:7,24,33,58,8:23,25:32,34:57,59:154)]

#convert from wide to long format
longtest1<-tidyr::gather(sorttest[ ,c(1:10,11:18)], "LUYR", "lu", c(11:18)) #this gathers properly, combine lu cols and retain identifiers
#if you rerun gather on the longtest1 it multiplies the number of rows, no easy option for mutiple gathers
#reshape is another option, but requires consistent and recognizable naming conventions
#http://stackoverflow.com/questions/25925556/gather-multiple-sets-of-columns-with-tidyr 

##I'll just do several indpependent gathers then join based on pointid and date, which I'll add before joining
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

longtest10<-tidyr::gather(sorttest[ ,c(2,83:90)], "ripcntyr", "rip_cnt", c(2:9)) #this gathers properly
longtest10<-mutate(longtest10, year= as.numeric(substr(longtest10$ripcntyr, 9,10))) #add a year column

longtest11<-tidyr::gather(sorttest[ ,c(2,91:98)], "precntyr", "pre_cnt", c(2:9)) #this gathers properly
longtest11<-mutate(longtest11, year= as.numeric(substr(longtest11$precntyr, 9,10))) #add a year column

longtest12<-tidyr::gather(sorttest[ ,c(2,99:106)], "agcntyr", "ag_cnt", c(2:9)) #this gathers properly
longtest12<-mutate(longtest12, year= as.numeric(substr(longtest12$agcntyr, 8,9))) #add a year column

longtest13<-tidyr::gather(sorttest[ ,c(2,107:114)], "domcntyr", "dom_cnt", c(2:9)) #this gathers properly
longtest13<-mutate(longtest13, year= as.numeric(substr(longtest13$domcntyr, 9,10))) #add a year column

longtest14<-tidyr::gather(sorttest[ ,c(2,115:122)], "indcntyr", "ind_cnt", c(2:9)) #this gathers properly
longtest14<-mutate(longtest14, year= as.numeric(substr(longtest14$indcntyr, 9,10))) #add a year column

longtest15<-tidyr::gather(sorttest[ ,c(2,123:130)], "fishcntyr", "fish_cnt", c(2:9)) #this gathers properly
longtest15<-mutate(longtest15, year= as.numeric(substr(longtest15$fishcntyr, 10,11))) #add a year column

longtest16<-tidyr::gather(sorttest[ ,c(2,131:138)], "reccntyr", "rec_cnt", c(2:9)) #this gathers properly
longtest16<-mutate(longtest16, year= as.numeric(substr(longtest16$reccntyr, 9,10))) #add a year column

longtest17<-tidyr::gather(sorttest[ ,c(2,139:146)], "adjudcntyr", "adjud_cnt", c(2:9)) #this gathers properly
longtest17<-mutate(longtest17, year= as.numeric(substr(longtest17$adjudcntyr, 11,12))) #add a year column

longtest18<-tidyr::gather(sorttest[ ,c(2,147:154)], "unauthcntyr", "unauth_cnt", c(2:9)) #this gathers properly
longtest18<-mutate(longtest18, year= as.numeric(substr(longtest18$unauthcntyr, 12,13))) #add a year column

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

#remove unneccessary columns
full.dat.long<-full.long[ ,c(-11,-14, -16,-18,-20,-22,-24,-26,-28,-30,-32,-34,-36,-38,-40,-42,-44,-46)]
#rearrange columns
full.dat.long<-full.dat.long[ ,c(1:10,12,11,13:29)]

####################################
#Add densities and percentages
####################################
#pixel res 1 km for calculating farmland area
full.dat.long <- mutate(full.dat.long, fl_areaKm = fl_count*1, wr_dens = wr_cnt/AreaSqKm, 
                        ag_dens = ag_cnt/fl_areaKm, rip_perc = rip_cnt/wr_cnt, pre_perc = pre_cnt/wr_cnt, 
                        ag_perc = ag_cnt/wr_cnt, dom_perc = dom_cnt/wr_cnt, ind_perc = ind_cnt/wr_cnt, 
                        fish_perc=fish_cnt/wr_cnt, rec_perc=rec_cnt/wr_cnt, adjud_perc = adjud_cnt/wr_cnt, unauth_perc=unauth_cnt/wr_cnt)

#############################
#Full long data as a Shapefile
###############################

cv.dat<-full.dat
cv.dat<- cv.dat[ ,2] #retain locations and POINTID only
cv.dat@data<-left_join(cv.dat@data,full.dat.long, by = "POINTID")

#plots to check joins
plot(cv.dat, cex = 0.005)

sub<-cv.dat
sub@data<-sub@data[sub@data$year==14,] #pull out one year at a time
spplot(sub[ ,"wr_cnt"], cex=0.005)
spplot(sub[ ,"wr_dens"], cex=0.005)
spplot(sub[ , "ag_cnt"], cex=0.005)
spplot(sub[ , "adjud_cnt"], cex=0.005)
spplot(sub[ , "gw"]) #looks good

#plot(cv.dat[cv.dat@data$year == 7 ,], cex =0.005) #indexing issue, proably becuase actually pull the complete 
#set of spatial locations for each eyar and it's expecting a subset

##########################################
#Write it out
##########################################
writeOGR(cv.dat, "/home/kate/CA/", "cv.dat", driver = "ESRI Shapefile")
saveRDS(cv.dat, file = "/home/kate/CA/data/cv.dat.rds")

#y<-readRDS('/home/kate/CA/data/cv.dat.rds') #check RDS save