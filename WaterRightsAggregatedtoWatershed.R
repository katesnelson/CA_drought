# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(maptools)

setwd("C:/")

#This script aggregates the cleaned water rights data created in WaterRightsDataClean.R to HUC12 watersheds. 
#To create a space-time dataset rerun the indented section for each year (changing the year label in object names)


##############################################
### IMPORT CENTRAL Valley HUC12 SHAPEFILE ####
##############################################

cv_huc <- readOGR(".","cv_huc")
cv_huc$huc12_id <- as.numeric(as.character(cv_huc$HUC12)) #fix the huc12_id which has duplicates
proj4string(cv_huc) #NAD83
cv_huc <- spTransform(cv_huc,CRS("+proj=longlat +ellps=WGS84"))
keep <-c("huc12_id","AreaAcres", "AreaSqKm", "diverse", "Name", "HUC12") #names of columns in cv_huc to keep
cv_huc <- cv_huc[, (names(cv_huc) %in% keep)] #reduce cv_huc to critical attributes


####Repeat this section for each year

    ######################################
    ### IMPORT WATER RIGHTS DATA  #######
    ######################################
    
    full.pt.2015<- readOGR(".","full.pt.2015") 
    full.pt.2015 <- spTransform(full.pt.2015,CRS("+proj=longlat +ellps=WGS84"))
    
    ##################################
    ###Associate WR Data with HUCs####
    ##################################
    
    pt.huc.2015 <- raster::intersect(full.pt.2015, cv_huc ) #spatial join of water right point of distribution (POD) pts and huc polygons
    pt.huc.dat <- pt.huc.2015@data #data only from join
    cv_huc_dup <-cv_huc #duplicate so we don't screw up the original
    
    ###########################################################
    ###Calculate aggregate values of interest for each HUC###
    ############################################################
    
    #set of placeholders
    wr_cnt<-0
    approp_cnt<-0
    junior_cnt<-0
    rip_cnt <-0
    pre_cnt<-0
    adjud_cnt<-0
    unauth_cnt<-0
    ag_cnt<-0
    dom_cnt<-0
    ind_cnt<-0
    fish_cnt<-0
    rec_cnt<-0
    
    #loop through hucs
    for (i in 1:nrow(cv_huc_dup)){
      subset<- dplyr::filter(pt.huc.dat, huc12_id == cv_huc_dup$huc12_id[i]) #select only wr PODs in a single HUC
      
      wr_cnt[i]<- length(subset$huc12_id)
      cv_huc_dup$wr_cnt[i] <- wr_cnt[i]  #add count of all water right PODs for each huc
      
      approp <- dplyr::filter(subset, WR.Type == "Appropriative")
      approp_cnt[i]<- length(approp$WR.Type)
      cv_huc_dup$approp_cnt[i] <- approp_cnt[i] #add count of appropriative wr PODs for each huc
      
      junior <- dplyr::filter(subset, is.na(Riparin) & is.na(Pre1914))
      junior_cnt[i]<- length(junior$WR.Type)
      cv_huc_dup$junior_cnt[i] <- junior_cnt[i] #add count of junior wr PODs for each huc
      
      rip <- dplyr::filter(subset, !is.na(Riparin))
      rip_cnt[i]<- length(rip$Riparin)
      cv_huc_dup$rip_cnt[i] <- rip_cnt[i] #add count of riparian wr PODs for each huc
      
      pre <- dplyr::filter(subset, !is.na(Pre1914))
      pre_cnt[i]<- length(pre$Pre1914)
      cv_huc_dup$pre_cnt[i] <- pre_cnt[i] #add count of pre1914 wr PODs for each huc
      
      ag <- dplyr::filter(subset, use == "Agriculture" )
      ag_cnt[i]<- length(ag$use)
      cv_huc_dup$ag_cnt[i] <- ag_cnt[i] #add count of agricultural wr PODs for each huc
      
      dom <- dplyr::filter(subset, use == "Domestic" )
      dom_cnt[i]<- length(dom$use)
      cv_huc_dup$dom_cnt[i] <- dom_cnt[i] #add count of domestic wr PODs for each huc
      
      ind <- dplyr::filter(subset, use == "Industrial" )
      ind_cnt[i]<- length(ind$use)
      cv_huc_dup$ind_cnt[i] <- ind_cnt[i] #add count of industrial wr PODs for each huc
      
      fish <- dplyr::filter(subset, use == "Fish & Wildlife" )
      fish_cnt[i]<- length(fish$use)
      cv_huc_dup$fish_cnt[i] <- fish_cnt[i] #add count of fish & wildlife wr PODs for each huc
      
      rec <- dplyr::filter(subset, use == "Recreation & Other" )
      rec_cnt[i]<- length(rec$use)
      cv_huc_dup$rec_cnt[i] <- rec_cnt[i] #add count of recreational wr PODs for each huc
      
      adjud <- dplyr::filter(subset, !is.na(CortDcr))
      adjud_cnt[i]<- length(adjud$CortDcr)
      cv_huc_dup$adjud_cnt[i] <- adjud_cnt[i] #add count of court ajudicated wr PODs for each huc
      
      unauth <- dplyr::filter(subset, !is.na(Un.Athr))
      unauth_cnt[i]<- length(unauth$Un.Athr)
      cv_huc_dup$unauth_cnt[i] <- unauth_cnt[i] #add count of  unauthorized wr PODs for each huc
    }
    
    
    cv_huc_dup <- spTransform(cv_huc_dup,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") ) #set to NAD83
    writeOGR(cv_huc_dup,".","finalwrbyhuc.2015", driver="ESRI Shapefile") #write the water right by huc shapefile

    
###########################################################################
###READ in and MERGE all of the Water Right POD Aggregated to HUC Info ###
############################################################################

wrbyhuc.2006<- readOGR(".","finalwrbyhuc.2006")
wrbyhuc.2007<- readOGR(".","finalwrbyhuc.2007")
wrbyhuc.2008<- readOGR(".","finalwrbyhuc.2008")
wrbyhuc.2009<- readOGR(".","finalwrbyhuc.2009")
wrbyhuc.2010<- readOGR(".","finalwrbyhuc.2010")
wrbyhuc.2011<- readOGR(".","finalwrbyhuc.2011")
wrbyhuc.2012<- readOGR(".","finalwrbyhuc.2012")
wrbyhuc.2013<- readOGR(".","finalwrbyhuc.2013")
wrbyhuc.2014<- readOGR(".","finalwrbyhuc.2014")
wrbyhuc.2015<- readOGR(".","finalwrbyhuc.2015")

full.wrbyhuc <- sp::merge(wrbyhuc.2006, wrbyhuc.2007@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2008@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2009@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2010@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2011@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2012@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2013@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2014@data, by = "huc12_id")
full.wrbyhuc <- sp::merge(full.wrbyhuc, wrbyhuc.2015@data, by = "huc12_id")

names(full.wrbyhuc) = c("huc12_id.06", "AreaAcres.06",  "AreaSqKm.06",   "HUC12_id",      "Name.06",    "diverse.06",    "wr_cnt.06", "approp_cnt.06", "junior_cnt.06",   "rip_cnt.06",    "pre_cnt.06",   
                       "ag_cnt.06",     "dom_cnt.06",    "ind_cnt.06",    "fish_cnt.06",   "rec_cnt.06",    "adjud_cnt.06",  "unauth_cnt.06",
                       "AreaAcres.07",  "AreaSqKm.07",    "HUC12_id.07",    "Name.07",       "diverse.07",    "wr_cnt.07",  "approp_cnt.07", "junior_cnt.07",   "rip_cnt.07",    "pre_cnt.07",   
                       "ag_cnt.07",     "dom_cnt.07",    "ind_cnt.07",    "fish_cnt.07",   "rec_cnt.07",    "adjud_cnt.07",  "unauth_cnt.07",
                       "AreaAcres.08",  "AreaSqKm.08",     "HUC12_id.08",    "Name.08",         "diverse.08",    "wr_cnt.08",  "approp_cnt.08", "junior_cnt.08",   "rip_cnt.08",    "pre_cnt.08",   
                       "ag_cnt.08",     "dom_cnt.08",    "ind_cnt.08",    "fish_cnt.08",   "rec_cnt.08",    "adjud_cnt.08",  "unauth_cnt.08",
                       "AreaAcres.09",  "AreaSqKm.09",     "HUC12_id.09", "Name.09",        "diverse.09",    "wr_cnt.09", "approp_cnt.09", "junior_cnt.09",    "rip_cnt.09",    "pre_cnt.09",   
                       "ag_cnt.09",     "dom_cnt.09",    "ind_cnt.09",    "fish_cnt.09",   "rec_cnt.09",    "adjud_cnt.09",  "unauth_cnt.09",
                       "AreaAcres.10",  "AreaSqKm.10",      "HUC12_id.10",  "Name.10",         "diverse.10",    "wr_cnt.10",  "approp_cnt.10", "junior_cnt.10",   "rip_cnt.10",    "pre_cnt.10",   
                       "ag_cnt.10",     "dom_cnt.10",    "ind_cnt.10",    "fish_cnt.10",   "rec_cnt.10",    "adjud_cnt.10",  "unauth_cnt.10",
                       "AreaAcres.11",  "AreaSqKm.11",      "HUC12_id.11",   "Name.11",        "diverse.11",    "wr_cnt.11",  "approp_cnt.11", "junior_cnt.11",   "rip_cnt.11",    "pre_cnt.11",   
                       "ag_cnt.11",     "dom_cnt.11",    "ind_cnt.11",    "fish_cnt.11",   "rec_cnt.11",    "adjud_cnt.11",  "unauth_cnt.11",
                       "AreaAcres.12",  "AreaSqKm.12",      "HUC12_id.12",  "Name.12",        "diverse.12",    "wr_cnt.12",   "approp_cnt.12", "junior_cnt.12",  "rip_cnt.12",    "pre_cnt.12",   
                       "ag_cnt.12",     "dom_cnt.12",    "ind_cnt.12",    "fish_cnt.12",   "rec_cnt.12",    "adjud_cnt.12",  "unauth_cnt.12",
                       "AreaAcres.13",  "AreaSqKm.13",      "HUC12_id.13", "Name.13",         "diverse.13",    "wr_cnt.13",    "approp_cnt.13", "junior_cnt.13", "rip_cnt.13",    "pre_cnt.13",   
                       "ag_cnt.13",     "dom_cnt.13",    "ind_cnt.13",    "fish_cnt.13",   "rec_cnt.13",    "adjud_cnt.13",  "unauth_cnt.13",
                       "AreaAcres.14",  "AreaSqKm.14",      "HUC12_id.14", "Name.14",          "diverse.14",    "wr_cnt.14",  "approp_cnt.14", "junior_cnt.14",   "rip_cnt.14",    "pre_cnt.14",   
                       "ag_cnt.14",     "dom_cnt.14",    "ind_cnt.14",    "fish_cnt.14",   "rec_cnt.14",    "adjud_cnt.14",  "unauth_cnt.14",
                       "AreaAcres.15",  "AreaSqKm.15",      "HUC12_id.15",  "Name.15",         "diverse.15",    "wr_cnt.15",   "approp_cnt.15", "junior_cnt.15",  "rip_cnt.15",    "pre_cnt.15",   
                       "ag_cnt.15",     "dom_cnt.15",    "ind_cnt.15",    "fish_cnt.15",   "rec_cnt.15",    "adjud_cnt.15",  "unauth_cnt.15")

#writeOGR(full.wrbyhuc,".","finalfull.wrbyhuc", driver="ESRI Shapefile") #write the farmland shapefile

saveRDS(full.wrbyhuc, file = "finalfullwrbyhuc.rds")
