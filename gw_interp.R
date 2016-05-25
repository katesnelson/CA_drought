### Game Plan ###

#Drop obs before 2000
#Yr (4 digit)-months(2 digit) as columns
#If multiple per months drop first (keep last)
#Wells as rows
#obs of elevation in the cells

#Then we want to interpolate to get value at each year
#Calculate delta from year to year
#Take average for hucs

#Could also do average for groundwater basins
#or mesh of points

###Start Script ###

# Required packages
library(dplyr)
library(sp)
library(rgdal)
library(spdep)
library(taRifx.geo)
library(rgeos)
library(raster)
library(readxl)

#setwd("C:/Users/nelsonks/Documents/Research/Soc Sci/Drought-Vul Data/central valley")
setwd("C:/Users/tuan/Documents/Research")
#setwd("C:/Users/tuan/Dropbox/My working papers")

#read in the gw point shapefile, tidy it up, and extract data to a  dataframe
gw_pts <- readOGR(".","gw_obs") 


names(gw_pts@data)<-c("SOURCE", "WELL.NUMBER", "MEASUREMENT.DATE", "DEPTH.TO.WATER", "GW.ELEVATION", "CASGEM.Well.Number","Local.Well.Designation", "Authority.Type",                                
                     "Co.operating.Agency", "Groundwater.Basin.Subbasin.Name","Groundwater.Basin..Subbasin.Number",            
                     "County", "Type.of.Well","Status.of.Well", "Well.Usage","Total.Well.Depth",  "Measurement.Count",                             
                     "Earliest.Elevation.Measurement.Date", "Most.Recent.Elevation.Measurement.Date",  
                     "Minimum.Groundwater.Elevation.Measured",         "Minimum.Groundwater.Elevation.Measurement.Date",
                     "Maximum.Groundwater.Elevation.Measured",         "Maximum.Groundwater.Elevation.Measurement.Date",
                     "Latitude..NAD.83.", "Longitude..NAD.83."  )
gw_pts@data = transform(gw_pts@data,
                   SOURCE= as.character(SOURCE),
                   WELL.NUMBER = as.character(WELL.NUMBER),
                   MEASUREMENT.DATE = as.character(MEASUREMENT.DATE),
                   DEPTH.TO.WATER = as.numeric(DEPTH.TO.WATER),
                   GW.ELEVATION = as.numeric(GW.ELEVATION),
                       CASGEM.Well.Number = as.character(CASGEM.Well.Number),
                       Local.Well.Designation = as.character(Local.Well.Designation),
                       Authority.Type = as.character(Authority.Type),
                       Co.operating.Agency = as.character(Co.operating.Agency),
                       Groundwater.Basin.Subbasin.Name = as.character(Groundwater.Basin.Subbasin.Name),
                       Groundwater.Basin..Subbasin.Number = as.character(Groundwater.Basin..Subbasin.Number),
                       County = as.character(County),
                       Type.of.Well = as.character(Type.of.Well),
                       Status.of.Well = as.character(Status.of.Well),
                       Well.Usage = as.character(Well.Usage),
                       Total.Well.Depth = as.numeric(Total.Well.Depth),
                       Measurement.Count = as.numeric(Measurement.Count),
                     Earliest.Elevation.Measurement.Date = as.character(Earliest.Elevation.Measurement.Date),
                       Most.Recent.Elevation.Measurement.Date= as.character(Most.Recent.Elevation.Measurement.Date),
                       Minimum.Groundwater.Elevation.Measured = as.numeric(Minimum.Groundwater.Elevation.Measured),
                       Maximum.Groundwater.Elevation.Measured = as.numeric(Maximum.Groundwater.Elevation.Measured),
                       Latitude..NAD.83.= as.numeric(Latitude..NAD.83.),
                       Longitude..NAD.83.= as.numeric(Longitude..NAD.83.)
)

gw_pts@data = transform(gw_pts@data,
                        MEASUREMENT.DATE = as.Date(MEASUREMENT.DATE),
                        Earliest.Elevation.Measurement.Date = as.Date(Earliest.Elevation.Measurement.Date),
                        Most.Recent.Elevation.Measurement.Date = as.Date(Most.Recent.Elevation.Measurement.Date)
)

gw_pts@data$yr.mnth <- format(gw_pts@data$MEASUREMENT.DATE, format = "%Y-%m") #create a column of 4-digit year and 2-digit month

#gw_dat <- gw_pts@data
#gw_dat$MEASUREMENT.DATE <- as.Date(gw_dat$MEASUREMENT.DATE)
#gw_dat$yr.mnth <- format(gw_dat$MEASUREMENT.DATE, format = "%Y-%m") #create a column of 4-digit year and 2-digit month


#create a new well observation time series dataframe
wells <- as.data.frame(gw_pts@data$WELL.NUMBER) #wells as rows
names(wells)<-c("WELL.NUMBER")
wells$WELL.NUMBER <- as.character(wells$WELL.NUMBER)
wells<-distinct(wells) # remove duplicate well records

dates <- seq(as.Date("2000-01-01"), by="month", length = 12*15) #create a sequence of dates
dates <- format(dates, format = "%Y-%m")                         #look only at year-month for dates
empty <- as.data.frame(matrix(data=NA, nrow = length(wells[,1]), ncol = length(dates))) #create an empty dataframe of length dates
colnames(empty) <- dates #use dates as column names of empty dataframe
ts <-bind_cols(wells, empty) #combine the dataframes to build the empty time series

#ts<- read.csv(".","gw_timeser.csv") 

for (i in 1:nrow(ts)){ #restart at ~7627, this loop is slow slow slow
  for (j in 1: length(dates)){
            well <- ts[i,1]
            date <- colnames(ts)[j]
            record <- filter(gw_pts@data, WELL.NUMBER %in% well)
            record <- filter(record, yr.mnth %in% date)
            if (nrow(record) != 0) {
           ts[i,j] <- first(record$GW.ELEVATION) 
  }
}
}

write.csv(ts,"gw_timeser.csv")



