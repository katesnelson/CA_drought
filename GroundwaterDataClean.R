# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(readxl)

setwd("C:/")

#This script reads in and creates a combined shapefile of California groundWater well data from csv files downloaded from the 
#California Statewide Groundwater Elevation Monitoring Program (CASGEM) at http://water.ca.gov/groundwater/casgem/index.cfm
#and the California State Water Resources Control Board Groundwater Ambient Monitoring and Assessment Program (GAMA) GeoTracker Online Database
#at http://geotracker.waterboards.ca.gov/gama/data_download.asp. 
#To use the script download the data, taking note of the file naming conventions used below.


#######################
#IMPORT CASGEM DATA ###
#######################
#Groundwater elevation data can be difficult to extract from CASGEM (the download process is slow and prone to disconnection).
#Therefore only well descriptive information, for both official CASGEM wells and voluntary wells, was downloaded.

CASGEM_wells <- read.csv("WellsCASGEM.csv")
CASGEM_vol <- read.csv("cagsem_wells_voluntary.csv")

################################
#IMPORT Geotracker GAMA DATA ###
################################
# The GAMA Data Download site includes a download link for Statewide Depth-to-Water and Groundwater Elevation data which is the most comprehensive and easy to 
#download source of California state-wide groundwater elevation measurements.According to a GeoTracker GAMA factsheet 
#(http://www.waterboards.ca.gov/publications_forms/publications/factsheets/docs/geotrkgama_fs_2015oct.pdf)this data 
#comes from Water Boards cleanup sites and the Department of Water Resources water data library (http://www.water.ca.gov/waterdatalibrary/groundwater/index.cfm).

gama <- read.csv("gama_all_dtw_elev.txt", sep="\t")

#Clean up GAMA and CASGEM data types and names

CASGEM_vol <- CASGEM_vol[1:22] #remove empty columns in the voluntary well descriptions object

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

#join CASGEM information with readings from GAMA by State Well Number
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


###Write Groundwater point data to ESRI shapefiles ###
writeOGR(gw_obs,".","gw_obs", driver="ESRI Shapefile")
