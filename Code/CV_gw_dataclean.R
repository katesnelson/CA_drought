# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(readxl)

setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/central valley")
#setwd("C:/Users/tuan/Documents/Research")
#setwd("C:/Users/tuan/Dropbox/My working papers")

############################################################################
### IMPORT Ground water data and clean and join to make more managable ### 
#####################################################################################

#####################################################
#GICIMA DATA, 
#https://gis.water.ca.gov/app/gicima/
#####################################################
#these files from GICIMA give gw points and elevation for wells with reading at both yrs 
#useful only for the early year density
#will not pick up wells that were active in earlier yr and dried up before the later yr
#NOT IN USE

#S2012_S2007 <- readOGR("S2012_S2007_Change_Points","S2012_05yr_pts") 
#S2013_S2008 <- readOGR("S013_S2008_Change_Points","S2013_05yr_pts")
#S2014_S2009 <- readOGR("S2014_S2009_Change_Points","S2014_05yr_pts")
#S2015_S2010 <- readOGR("S2015_S2010_Change_Points","S2015_05yr_pts")
#S2014_S2011 <- readOGR("S2014_S2011_Change_Points","S2014_03yr_pts")
#S2015_S2012 <- readOGR("S2015_S2012_Change_Points","S2015_03yr_pts")
#S2014_S2013 <- readOGR("S2014_S2013_Change_Points","S2014_01yr_pts")
#S2015_S2014 <- readOGR("S2015_S2014_Change_Points","S2015_01yr_pts")
#S2015 <- readOGR("S2015_DGBS_Points")

#dat<-as.data.frame(S2012_S2007@data)
#names07 <- c("ObjectID","SiteID","Depth07","Depth12","Diff","EWSEL","LWSEL","Date07","EWLM","Date12","LWLM")
#names(S2012_S2007) <- names07
#keep07 <- c("SiteID","Depth07","Date07")
#GW2007 <- S2012_S2007[, (names(S2012_S2007)%in% keep07)]


############################################################################
#CASGEM California Statewide Groundwater Elevation Monitoring Program
#http://water.ca.gov/groundwater/casgem/index.cfm
#############################################################################
#data from CASGEM, readings go until 9/9/2015, 
#this does not provide individual gw elevation readings, but provides well information (including location)
#gw elevation reading can be obtained from CASGEM, but the download process is slow and prone to disconnection

CASGEM_wells <- read.csv("WellsCASGEM.csv")
max(as.character(CASGEM_wells$Most.Recent.Elevation.Measurement.Date))
CASGEM_vol <- read.csv("cagsem_wells_voluntary.csv")
max(as.character(CASGEM_vol$Most.Recent.Elevation.Measurement.Date))

##############################################################################
#Geotracker GAMA 
#http://geotracker.waterboards.ca.gov/gama/data_download.asp
##################################################################################
# (as enar as i can tell this elev data is a combination of CASGEM, USGS, and GAMA observations), readings until 12/31/2014
# http://www.waterboards.ca.gov/publications_forms/publications/factsheets/docs/geotrkgama_fs_2015oct.pdf
gama <- read.csv("gama_all_dtw_elev.txt", sep="\t")
max(as.character(gama$MEASUREMENT.DATE))


#Clean up gama and casgem data types and names

CASGEM_vol <- CASGEM_vol[1:22] #remove empty columns

names(CASGEM_vol)<-c("State.Well.Number", "CASGEM.Well.Number","Local.Well.Designation", "Authority.Type",                                
                     "Monitoring.Entity", "Co.operating.Agency", "Groundwater.Basin.Subbasin.Name","Groundwater.Basin..Subbasin.Number",            
                     "County", "Type.of.Well","Status.of.Well", "Well.Usage","Total.Well.Depth",  "Measurement.Count",                             
                     "Earliest.Elevation.Measurement.Date", "Most.Recent.Elevation.Measurement.Date",  
                     "Minimum.Groundwater.Elevation.Measured",         "Minimum.Groundwater.Elevation.Measurement.Date",
                     "Maximum.Groundwater.Elevation.Measured",         "Maximum.Groundwater.Elevation.Measurement.Date",
                     "Latitude..NAD.83.", "Longitude..NAD.83."  )
CASGEM_vol = transform(CASGEM_vol, 
                       State.Well.Number = as.character(State.Well.Number),
                       CASGEM.Well.Number = as.character(CASGEM.Well.Number),
                       Local.Well.Designation = as.character(Local.Well.Designation),
                       Authority.Type = as.character(Authority.Type),
                       Monitoring.Entity = as.character(Monitoring.Entity),
                       Co.operating.Agency = as.character(Co.operating.Agency),
                       Groundwater.Basin.Subbasin.Name = as.character(Groundwater.Basin.Subbasin.Name),
                       Groundwater.Basin..Subbasin.Number = as.character(Groundwater.Basin..Subbasin.Number),
                       County = as.character(County),
                       Type.of.Well = as.character(Type.of.Well),
                       Status.of.Well = as.character(Status.of.Well),
                       Well.Usage = as.character(Well.Usage),
                       Total.Well.Depth = as.numeric(Total.Well.Depth),
                       Measurement.Count = as.numeric(Measurement.Count),
                       Earliest.Elevation.Measurement.Date = as.Date(Earliest.Elevation.Measurement.Date, "%m/%d/%Y"),
                       Most.Recent.Elevation.Measurement.Date= as.Date(Most.Recent.Elevation.Measurement.Date, "%m/%d/%Y"),
                       Minimum.Groundwater.Elevation.Measured = as.numeric(Minimum.Groundwater.Elevation.Measured),
                       Minimum.Groundwater.Elevation.Measurement.Date = as.Date(Minimum.Groundwater.Elevation.Measurement.Date, "%m/%d/%Y"),
                       Maximum.Groundwater.Elevation.Measured = as.numeric(Maximum.Groundwater.Elevation.Measured),
                       Maximum.Groundwater.Elevation.Measurement.Date = as.Date(Maximum.Groundwater.Elevation.Measurement.Date, "%m/%d/%Y"),
                       Latitude..NAD.83.= as.numeric(Latitude..NAD.83.),
                       Longitude..NAD.83.= as.numeric(Longitude..NAD.83.)
)


names(CASGEM_wells)<-c("State.Well.Number", "CASGEM.Well.Number","Local.Well.Designation", "Authority.Type",                                
                       "Monitoring.Entity", "Co.operating.Agency", "Groundwater.Basin.Subbasin.Name","Groundwater.Basin..Subbasin.Number",            
                       "County", "Type.of.Well","Status.of.Well", "Well.Usage","Total.Well.Depth",  "Measurement.Count",                             
                       "Earliest.Elevation.Measurement.Date", "Most.Recent.Elevation.Measurement.Date",  
                       "Minimum.Groundwater.Elevation.Measured",         "Minimum.Groundwater.Elevation.Measurement.Date",
                       "Maximum.Groundwater.Elevation.Measured",         "Maximum.Groundwater.Elevation.Measurement.Date",
                       "Latitude..NAD.83.", "Longitude..NAD.83."  )
CASGEM_wells = transform(CASGEM_wells, 
                         State.Well.Number = as.character(State.Well.Number),
                         CASGEM.Well.Number = as.character(CASGEM.Well.Number),
                         Local.Well.Designation = as.character(Local.Well.Designation),
                         Authority.Type = as.character(Authority.Type),
                         Monitoring.Entity = as.character(Monitoring.Entity),
                         Co.operating.Agency = as.character(Co.operating.Agency),
                         Groundwater.Basin.Subbasin.Name = as.character(Groundwater.Basin.Subbasin.Name),
                         Groundwater.Basin..Subbasin.Number = as.character(Groundwater.Basin..Subbasin.Number),
                         County = as.character(County),
                         Type.of.Well = as.character(Type.of.Well),
                         Status.of.Well = as.character(Status.of.Well),
                         Well.Usage = as.character(Well.Usage),
                         Total.Well.Depth = as.numeric(Total.Well.Depth),
                         Measurement.Count = as.numeric(Measurement.Count),
                         Earliest.Elevation.Measurement.Date = as.Date(Earliest.Elevation.Measurement.Date, "%m/%d/%Y"),
                         Most.Recent.Elevation.Measurement.Date= as.Date(Most.Recent.Elevation.Measurement.Date, "%m/%d/%Y"),
                         Minimum.Groundwater.Elevation.Measured = as.numeric(Minimum.Groundwater.Elevation.Measured),
                         Minimum.Groundwater.Elevation.Measurement.Date = as.Date(Minimum.Groundwater.Elevation.Measurement.Date, "%m/%d/%Y"),
                         Maximum.Groundwater.Elevation.Measured = as.numeric(Maximum.Groundwater.Elevation.Measured),
                         Maximum.Groundwater.Elevation.Measurement.Date = as.Date(Maximum.Groundwater.Elevation.Measurement.Date, "%m/%d/%Y"),
                         Latitude..NAD.83.= as.numeric(Latitude..NAD.83.),
                         Longitude..NAD.83.= as.numeric(Longitude..NAD.83.)
)

gama <- transform(gama,
                  SOURCE= as.character(SOURCE),
                  WELL.NUMBER = as.character(WELL.NUMBER),
                  MEASUREMENT.DATE = as.Date(MEASUREMENT.DATE, "%m/%d/%Y"),
                  DEPTH.TO.WATER = as.numeric(DEPTH.TO.WATER),
                  GW.ELEVATION = as.numeric(GW.ELEVATION),
                  LATITUDE = as.numeric(LATITUDE),
                  LONGITUDE = as.numeric(LONGITUDE)
)

#merge the official and voluntary CASGEM well records
CASGEM <- bind_rows (CASGEM_vol, CASGEM_wells)

#join CASGEM information with readings from gama by State Well Number
gama_CASGEM <-left_join(gama, CASGEM, by = c("WELL.NUMBER" = "State.Well.Number"))


#drop obs before 2000 and duplicate rows
gw_obs <- gama_CASGEM[!(gama_CASGEM$'MEASUREMENT.DATE' < as.Date("2000-01-01")),] 
gw_obs <- gw_obs[(gw_obs$'SOURCE' == 'DWR Well' ),] 
gw_obs <- gw_obs[!is.na(gw_obs$'WELL.NUMBER' ),] 
gw_obs <- distinct(gw_obs)
dup_chk <- gw_obs[-11] #remove monitoring entity column as it produces duplicates for joint monitored wells
dup_chk <- distinct(dup_chk)
gw_obs <- dup_chk

###Convert to SpatialPointsDataFrame and Remove Spatial Duplicates###
coordinates(gw_obs) = c('LONGITUDE','LATITUDE')#set coordinate system
proj4string(gw_obs)<- "+proj=longlat +ellps=WGS84"
plot(gw_obs, pch = 16, col = "red", cex = 0.3)



###Write gw point data to ESRI shapefiles ###
writeOGR(gw_obs,".","gw_obs", driver="ESRI Shapefile")
