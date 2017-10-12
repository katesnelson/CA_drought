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

#This script reads in California Water Rights Data from excel files extracted from the California State Water Resources Control Board
#electronic water rights information management system (eWRIMS) accessed at http://www.waterboards.ca.gov/waterrights/water_issues/programs/ewrims/index.shtml
#and creates a shapefile for active water rights in California for a series of specified years. 
#Due to download size constraints data for the entire state was extracted by individual use category.
#To use the script download the data from the eWRIMS system in excel format, taking note of the file naming convention used below.
#Some data sheets produced errors when imported as excel files. These sheets (POD & USES sheets) should be exported 
#to a csv file from within excel prior to running the script.

################################################################
### IMPORT WATER RIGHTS DATA SHEETS FOR EACH USE TYPE & JOIN ### 
################################################################

##read in all irrigation water rights data
irr_wr <-read_excel("irrigation.xls", sheet = "Water Rights", col_names=TRUE) #read in the water rights sheet
irr_pod <-read.csv("irr_pod.csv") #read in point of distribution sheet #some problems with this sheet, export to a csv and read in seperately
irr_app <- read_excel("irrigation.xls", sheet = "Application Info", col_names=TRUE, #read in application info sheet, specify data type
           col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
irr_use <- read.csv("irr_uses.csv") #read in water uses sheet #some problems with this sheet, export to a csv and read in seperately

    #join all irrigation water rights information 
    irr_wr <- left_join (irr_wr, irr_use, by = c("Application Number" = "Application.Number")) #join use info to wr info by number
    irr_wrapp <- left_join(irr_wr, irr_app, by = c("Application Number" = "Application ID")) #join app info to wr & use info
    irr_podwrapp <- left_join(irr_pod, irr_wrapp, by = c("Application.Number" = "Application Number")) #join point of distribution info to other wr info
    write.csv(irr_podwrapp, "irr_podwrapp.csv") #write the full dataset for each use to its own csv
    rm(irr_wr, irr_pod, irr_app, irr_use, irr_wrapp) #clear out unnecessary info

##read in all domestic water rights data
dom_wr <-read_excel("domestic.xls", sheet = "Water Rights", col_names=TRUE)
dom_pod <-read.csv("dom_pod.csv") #some problems with this sheet, export to a csv and read in seperately
dom_app <- read_excel("domestic.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
dom_use <- read.csv("dom_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join all domestic water rights information
    dom_wr <- left_join (dom_wr, dom_use, by = c("Application Number" = "Application.Number") )
    dom_wrapp <- left_join(dom_wr, dom_app, by = c("Application Number" = "Application ID"))
    dom_podwrapp <- left_join(dom_pod, dom_wrapp, by = c("Application.Number" = "Application Number"))
    write.csv(dom_podwrapp, "dom_podwrapp.csv")
    rm(dom_wr, dom_pod, dom_app, dom_use, dom_wrapp) #clear out unnecessary info

##read in all industrial water right data
ind_wr <-read_excel("industrial.xls", sheet = "Water Rights", col_names=TRUE)
ind_pod <-read.csv("ind_pod.csv")
ind_app <- read_excel("industrial.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
ind_use <- read.csv("ind_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join all industrial water right information
    ind_wr <- left_join (ind_wr, ind_use, by = c("Application Number" = "Application.Number") )
    ind_wrapp <- left_join(ind_wr, ind_app, by = c("Application Number" = "Application ID"))
    ind_podwrapp <- left_join(ind_pod, ind_wrapp, by = c("Application.Number" = "Application Number"))
    write.csv(ind_podwrapp, "ind_podwrapp.csv")
    rm(ind_wr, ind_pod, ind_app, ind_use, ind_wrapp) #clear out unnecessary info

##read in full stockwatering data
stk_wr <-read_excel("stockwater.xls", sheet = "Water Rights", col_names=TRUE)
stk_pod <-read.csv("stk_pod.csv")
stk_app <- read_excel("stockwater.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
stk_use <- read.csv("stk_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join stockwatering information
    stk_wr <- left_join (stk_wr, stk_use, by = c("Application Number" = "Application.Number") )
    stk_wrapp <- left_join(stk_wr, stk_app, by = c("Application Number" = "Application ID"))
    stk_podwrapp <- left_join(stk_pod, stk_wrapp, by = c("Application.Number" = "Application Number"))
    write.csv(stk_podwrapp, "stk_podwrapp.csv")
    rm(stk_wr, stk_pod, stk_app, stk_use, stk_wrapp) #clear out unnecessary info

##read in full power production data
pow_wr <-read_excel("power.xls", sheet = "Water Rights", col_names=TRUE)
pow_pod <-read.csv("pow_pod.csv")
pow_app <- read_excel("power.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
pow_use <- read.csv("pow_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join power production information
    pow_wr <- left_join (pow_wr, pow_use, by = c("Application Number" = "Application.Number") )
    pow_wrapp <- left_join(pow_wr, pow_app, by = c("Application Number" = "Application ID"))
    pow_podwrapp <- left_join(pow_pod, pow_wrapp, by = c("Application.Number" = "Application Number"))
    pow_podwrapp$POD.Num <- as.character(pow_podwrapp$POD.Num) #fix misspecified data classes
    pow_podwrapp$`Permit ID` <- as.character(pow_podwrapp$`Permit ID`)
    pow_podwrapp$DD.Beg.Month..Day <- as.character(pow_podwrapp$DD.Beg.Month..Day)
    pow_podwrapp$DD.End.Month.Day <- as.character(pow_podwrapp$DD.End.Month.Day)
    pow_podwrapp$Store.Beg.Month.Day <- as.character(pow_podwrapp$Store.Beg.Month.Day)
    write.csv(pow_podwrapp, "pow_podwrapp.csv")
    rm(pow_wr, pow_pod, pow_app, pow_use, pow_wrapp) #clear out unnecessary info

## read in full fire fighting water right data
fire_wr <-read_excel("fire.xls", sheet = "Water Rights", col_names=TRUE)
fire_pod <-read.csv("fire_pod.csv")
fire_app <- read_excel("fire.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
fire_use <- read.csv("fire_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join fire fighting information 
    fire_wr <- left_join (fire_wr, fire_use, by = c("Application Number" = "Application.Number") )
    fire_wrapp <- left_join(fire_wr, fire_app, by = c("Application Number" = "Application ID"))
    fire_podwrapp <- left_join(fire_pod, fire_wrapp, by = c("Application.Number" = "Application Number"))
    write.csv(fire_podwrapp, "fire_podwrapp.csv")
    rm(fire_wr, fire_pod, fire_app, fire_use, fire_wrapp) #clear out unnecessary info

##read in full wildlife protection water right data
wild_wr <-read_excel("wildlife.xls", sheet = "Water Rights", col_names=TRUE)
wild_pod <-read.csv("wild_pod.csv")
wild_app <- read_excel("wildlife.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
wild_use <- read.csv("wild_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join wildlife protection information
    wild_wr <- left_join (wild_wr, wild_use, by = c("Application Number" = "Application.Number") )
    wild_wrapp <- left_join(wild_wr, wild_app, by = c("Application Number" = "Application ID"))
    wild_podwrapp <- left_join(wild_pod, wild_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(wild_podwrapp, "wild_podwrapp.csv")
    rm(wild_wr, wild_pod, wild_app, wild_use, wild_wrapp) #clear out unnecessary info

##read in full frost protection water right data
fro_wr <-read_excel("frost.xls", sheet = "Water Rights", col_names=TRUE)
fro_pod <-read.csv("fro_pod.csv")
fro_app <- read_excel("frost.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
fro_use <- read.csv("fro_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join frost protection information
    fro_wr <- left_join (fro_wr, fro_use, by = c("Application Number" = "Application.Number") )
    fro_wrapp <- left_join(fro_wr, fro_app, by = c("Application Number" = "Application ID"))
    fro_podwrapp <- left_join(fro_pod, fro_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(fro_podwrapp, "fro_podwrapp.csv")
    rm(fro_wr, fro_pod, fro_app, fro_use, fro_wrapp) #clear out unnecessary info

##read in full heat protection water right data
heat_wr <-read_excel("heat.xls", sheet = "Water Rights", col_names=TRUE)
heat_pod <-read.csv("heat_pod.csv")
heat_app <- read_excel("heat.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
heat_use <- read.csv("heat_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join heat protection information
    heat_wr <- left_join (heat_wr, heat_use, by = c("Application Number" = "Application.Number") )
    heat_wrapp <- left_join(heat_wr, heat_app, by = c("Application Number" = "Application ID"))
    heat_podwrapp <- left_join(heat_pod, heat_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(heat_podwrapp, "heat_podwrapp.csv")
    rm(heat_wr, heat_pod, heat_app, heat_use, heat_wrapp) #clear out unnecessary info

##read in full incidental power water right data
incp_wr <-read_excel("incidental power.xls", sheet = "Water Rights", col_names=TRUE)
incp_pod <-read.csv("incp_pod.csv")
incp_app <- read_excel("incidental power.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
incp_use <- read.csv("incp_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join indcidental power information
    incp_wr <- left_join (incp_wr, incp_use, by = c("Application Number" = "Application.Number") )
    incp_wrapp <- left_join(incp_wr, incp_app, by = c("Application Number" = "Application ID"))
    incp_podwrapp <- left_join(incp_pod, incp_wrapp, by = c("Application.Number"="Application Number"))
    incp_podwrapp$POD.Num <- as.character(incp_podwrapp$POD.Num) #fix misspecified data classes
    incp_podwrapp$`Permit ID` <- as.character(incp_podwrapp$`Permit ID`)
    incp_podwrapp$DD.Beg.Month..Day <- as.character(incp_podwrapp$DD.Beg.Month..Day)
    incp_podwrapp$DD.End.Month.Day <- as.character(incp_podwrapp$DD.End.Month.Day)
    incp_podwrapp$Store.Beg.Month.Day <- as.character(incp_podwrapp$Store.Beg.Month.Day)
    write.csv(incp_podwrapp, "incp_podwrapp.csv")
    rm(incp_wr, incp_pod, incp_app, incp_use, incp_wrapp) #clear out unnecessary info

##read in full milling water rights data
mil_wr <-read_excel("milling.xls", sheet = "Water Rights", col_names=TRUE)
mil_pod <-read.csv("mil_pod.csv")
mil_app <- read_excel("milling.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
mil_use <- read.csv("mil_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join milling information
    mil_wr <- left_join (mil_wr, mil_use, by = c("Application Number" = "Application.Number") )
    mil_wrapp <- left_join(mil_wr, mil_app, by = c("Application Number" = "Application ID"))
    mil_podwrapp <- left_join(mil_pod, mil_wrapp, by = c("Application.Number"="Application Number"))
    mil_podwrapp$POD.Num <- as.character(mil_podwrapp$POD.Num) #fix misspecified data classes
    mil_podwrapp$`Permit ID` <- as.character(mil_podwrapp$`Permit ID`)
    mil_podwrapp$DD.Beg.Month..Day <- as.character(mil_podwrapp$DD.Beg.Month..Day)
    mil_podwrapp$DD.End.Month.Day <- as.character(mil_podwrapp$DD.End.Month.Day)
    mil_podwrapp$Store.Beg.Month.Day <- as.character(mil_podwrapp$Store.Beg.Month.Day)
    write.csv(mil_podwrapp, "mil_podwrapp.csv")
    rm(mil_wr, mil_pod, mil_app, mil_use, mil_wrapp) #clear out unnecessary info

##read in full mining water right data
min_wr <-read_excel("mining.xls", sheet = "Water Rights", col_names=TRUE)
min_pod <-read.csv("min_pod.csv")
min_app <- read_excel("mining.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
min_use <- read.csv("min_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join mining information
    min_wr <- left_join (min_wr, min_use, by = c("Application Number" = "Application.Number") )
    min_wrapp <- left_join(min_wr, min_app, by = c("Application Number" = "Application ID"))
    min_podwrapp <- left_join(min_pod, min_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(min_podwrapp, "min_podwrapp.csv")
    rm(min_wr, min_pod, min_app, min_use, min_wrapp) #clear out unnecessary info

##read in full municipal water right data
mun_wr <-read_excel("municipal.xls", sheet = "Water Rights", col_names=TRUE)
mun_pod <-read.csv("mun_pod.csv")
mun_app <- read_excel("municipal.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
mun_use <- read.csv("mun_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join municipal information
    mun_wr <- left_join (mun_wr, mun_use, by = c("Application Number" = "Application.Number") )
    mun_wrapp <- left_join(mun_wr, mun_app, by = c("Application Number" = "Application ID"))
    mun_podwrapp <- left_join(mun_pod, mun_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(mun_podwrapp, "mun_podwrapp.csv")
    rm(mun_wr, mun_pod, mun_app, mun_use, mun_wrapp) #clear out unnecessary info

##read in full other water right data
oth_wr <-read_excel("other.xls", sheet = "Water Rights", col_names=TRUE)
oth_pod <-read.csv("oth_pod.csv")
oth_app <- read_excel("other.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
oth_use <- read.csv("oth_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join other information 
    oth_wr <- left_join (oth_wr, oth_use, by = c("Application Number" = "Application.Number") )
    oth_wrapp <- left_join(oth_wr, oth_app, by = c("Application Number" = "Application ID"))
    oth_podwrapp <- left_join(oth_pod, oth_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(oth_podwrapp, "oth_podwrapp.csv")
    rm(oth_wr, oth_pod, oth_app, oth_use, oth_wrapp) #clear out unnecessary info

##read in full recreational water right data
rec_wr <-read_excel("recreation.xls", sheet = "Water Rights", col_names=TRUE)
rec_pod <-read.csv("rec_pod.csv")
rec_app <- read_excel("recreation.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
rec_use <- read.csv("rec_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join recreational information
    rec_wr <- left_join (rec_wr, rec_use, by = c("Application Number" = "Application.Number") )
    rec_wrapp <- left_join(rec_wr, rec_app, by = c("Application Number" = "Application ID"))
    rec_podwrapp <- left_join(rec_pod, rec_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(rec_podwrapp, "rec_podwrapp.csv")
    rm(rec_wr, rec_pod, rec_app, rec_use, rec_wrapp) #clear out unnecessary info

##read in full snow making water right data
snow_wr <-read_excel("snow making.xls", sheet = "Water Rights", col_names=TRUE)
snow_pod <-read.csv("snow_pod.csv")
snow_app <- read_excel("snow making.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
snow_use <- read.csv("snow_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join snow making information
    snow_wr <- left_join (snow_wr, snow_use, by = c("Application Number" = "Application.Number") )
    snow_wrapp <- left_join(snow_wr, snow_app, by = c("Application Number" = "Application ID"))
    snow_podwrapp <- left_join(snow_pod, snow_wrapp, by = c("Application.Number"="Application Number"))
    snow_podwrapp$POD.Num <- as.character(snow_podwrapp$POD.Num) #fix misspecifed data classes
    snow_podwrapp$`Permit ID` <- as.character(snow_podwrapp$`Permit ID`)
    snow_podwrapp$DD.Beg.Month..Day <- as.character(snow_podwrapp$DD.Beg.Month..Day)
    snow_podwrapp$DD.End.Month.Day <- as.character(snow_podwrapp$DD.End.Month.Day)
    snow_podwrapp$Store.Beg.Month.Day <- as.character(snow_podwrapp$Store.Beg.Month.Day)
    write.csv(snow_podwrapp, "snow_podwrapp.csv")
    rm(snow_wr, snow_pod, snow_app, snow_use, snow_wrapp) #clear out unnecessary info

##read in full aquaculuture water rights data
aqua_wr <-read_excel("aquaculture.xls", sheet = "Water Rights", col_names=TRUE)
aqua_pod <-read.csv("aqua_pod.csv")
aqua_app <- read_excel("aquaculture.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
aqua_use <- read.csv("aqua_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join aquaculture information
    aqua_wr <- left_join (aqua_wr, aqua_use, by = c("Application Number" = "Application.Number") )
    aqua_wrapp <- left_join(aqua_wr, aqua_app, by = c("Application Number" = "Application ID"))
    aqua_podwrapp <- left_join(aqua_pod, aqua_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(aqua_podwrapp, "aqua_podwrapp.csv")
    rm(aqua_wr, aqua_pod, aqua_app, aqua_use, aqua_wrapp) #clear out unnecessary info

##read in full aesthetic water right data
aes_wr <-read_excel("aesthetic.xls", sheet = "Water Rights", col_names=TRUE)
aes_pod <-read.csv("aes_pod.csv")
aes_app <- read_excel("aesthetic.xls", sheet = "Application Info", col_names=TRUE, 
                      col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
aes_use <- read.csv("aes_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    #join aesthetic information
    aes_wr <- left_join (aes_wr, aes_use, by = c("Application Number" = "Application.Number") )
    aes_wrapp <- left_join(aes_wr, aes_app, by = c("Application Number" = "Application ID"))
    aes_podwrapp <- left_join(aes_pod, aes_wrapp, by = c("Application.Number"="Application Number"))
    aes_podwrapp$POD.Num <- as.character(aes_podwrapp$POD.Num) #fix misspecified data classes
    aes_podwrapp$`Permit ID` <- as.character(aes_podwrapp$`Permit ID`)
    aes_podwrapp$DD.Beg.Month..Day <- as.character(aes_podwrapp$DD.Beg.Month..Day)
    aes_podwrapp$DD.End.Month.Day <- as.character(aes_podwrapp$DD.End.Month.Day)
    aes_podwrapp$Store.Beg.Month.Day <- as.character(aes_podwrapp$Store.Beg.Month.Day)
    write.csv(aes_podwrapp, "aes_podwrapp.csv")
    rm(aes_wr, aes_pod, aes_app, aes_use, aes_wrapp) #clear out unnecessary info

##read in full dust control water right data
dust_wr <-read_excel("dust.xls", sheet = "Water Rights", col_names=TRUE)
dust_pod <-read.csv("dust_pod.csv")
dust_app <- read_excel("dust.xls", sheet = "Application Info", col_names=TRUE, 
                       col_types = c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text"))
dust_use <- read.csv("dust_uses.csv") #some problems with this sheet, export to a csv and read in seperately

    # join dust control information
    dust_wr <- left_join (dust_wr, dust_use, by = c("Application Number" = "Application.Number") )
    dust_wrapp <- left_join(dust_wr, dust_app, by = c("Application Number" = "Application ID"))
    dust_podwrapp <- left_join(dust_pod, dust_wrapp, by = c("Application.Number"="Application Number"))
    write.csv(dust_podwrapp, "dust_podwrapp.csv")
    rm(dust_wr, dust_pod, dust_app, dust_use, dust_wrapp) #clear out unnecessary info

####################################################################
###JOIN INDIVIDUAL USE TYPE RECORDS IN SIMPLIFIED USE CATEGORIES####
####################################################################
    
#Conduct joins in order of relative importance (subjective decision) so that water rights with multiple uses
#are categorized under use of greatest importance 

#agriculture uses    
    ag <- full_join(irr_podwrapp,stk_podwrapp)
    ag <- full_join(ag,heat_podwrapp)
    ag <- full_join(ag,fro_podwrapp)
    ag <- full_join(ag,aqua_podwrapp)
    write.csv(ag, "agriculture_wr.csv")

#domestic uses
    dom <- full_join(dom_podwrapp,mun_podwrapp)
    dom <- full_join(dom,aes_podwrapp)
    write.csv(dom, "domestic_wr.csv") 

#industrial uses
    ind <- full_join(ind_podwrapp,dust_podwrapp)
    ind <- full_join(ind,pow_podwrapp)
    ind <- full_join(ind,min_podwrapp)
    ind <- full_join(ind,mil_podwrapp)
    ind <- full_join(ind,incp_podwrapp)
    write.csv(ind, "industrial_wr.csv") 

#fish & wildlife uses
    wild <- full_join(wild_podwrapp, fire_podwrapp)
    write.csv(wild, "wildlife_wr.csv") 

#recreational uses
    other <- full_join(rec_podwrapp,oth_podwrapp)
    other <- full_join(other,snow_podwrapp)
    write.csv(other, "other_wr.csv")  

#remove unneeded objects
rm(irr_podwrapp,stk_podwrapp,aqua_podwrapp,fro_podwrapp,heat_podwrapp,dom_podwrapp,mun_podwrapp,
   aes_podwrapp,ind_podwrapp, dust_podwrapp, pow_podwrapp,incp_podwrapp,min_podwrapp,mil_podwrapp,
   wild_podwrapp,fire_podwrapp,oth_podwrapp,rec_podwrapp,snow_podwrapp)

#######################################
### TRIM DATA TO MOST RELEVANT INFO ###
#######################################

keep <- c("Application.Number", "POD.ID", "Latitude", "Longitude", "POD.Status", "POD.Direct.Diversion.Rate", "DD.Unit", 
          "POD.Storage", "County.x", "POD.HUC.12.Name", "WR Type", "Status", "Status Date", "Face Value Amount", "Beneficial.Use",
          "Riparian", "Pre 1914", "Court Decree", "Un-Authorized")

ag_dat <- ag[,(names(ag) %in% keep)]    #trim to retain only the named columns
ag_dat$Latitude <- as.numeric(ag_dat$Latitude, options(digits = 15))    #format data for mapping
ag_dat$Longitude <- as.numeric(ag_dat$Longitude, options(digits = 15))    #format data for mapping
ag_dat$'Status Date' <- as.Date(ag_dat$'Status Date', "%m/%d/%Y")      #format dates
ag_dat <- filter(ag_dat, !is.na(Latitude))              #remove locations without spatial identifier information
ag_dat <- distinct(ag_dat) #removes duplicate records

dom_dat <- dom[,(names(dom) %in% keep)]
dom_dat$Latitude <- as.numeric(dom_dat$Latitude, options(digits = 15))
dom_dat$Longitude <- as.numeric(dom_dat$Longitude, options(digits = 15))
dom_dat$'Status Date' <- as.Date(dom_dat$'Status Date', "%m/%d/%Y")
dom_dat <- filter(dom_dat, !is.na(Latitude))
dom_dat <- distinct(dom_dat)

ind_dat <- ind[,(names(ind) %in% keep)]
ind_dat$Latitude <- as.numeric(ind_dat$Latitude, options(digits = 15))
ind_dat$Longitude <- as.numeric(ind_dat$Longitude, options(digits = 15))
ind_dat$'Status Date' <- as.Date(ind_dat$'Status Date', "%m/%d/%Y")
ind_dat <- filter(ind_dat, !is.na(Latitude))
ind_dat <- distinct(ind_dat)

wild_dat <- wild[,(names(wild) %in% keep)]
wild_dat$Latitude <- as.numeric(wild_dat$Latitude, options(digits = 15))
wild_dat$Longitude <- as.numeric(wild_dat$Longitude, options(digits = 15))
wild_dat$'Status Date' <- as.Date(wild_dat$'Status Date', "%m/%d/%Y")
wild_dat <- filter(wild_dat, !is.na(Latitude))
wild_dat <- distinct(wild_dat)

other_dat <- other[,(names(other) %in% keep)]
other_dat$Latitude <- as.numeric(other_dat$Latitude, options(digits = 15))
other_dat$Longitude <- as.numeric(other_dat$Longitude, options(digits = 15))
other_dat$'Status Date' <- as.Date(other_dat$'Status Date', "%m/%d/%Y")
other_dat <- filter(other_dat, !is.na(Latitude))
other_dat <- distinct(other_dat)

#####################################################################
###Convert to SpatialPointsDataFrame and Remove Spatial Duplicates###
######################################################################

#This process removes duplicate records for a single water right
#(eg. Point of diversion with multiple records associated with supplemental 
#statements of diversion required by CA SWRCB every 3 years 
#http://www.waterboards.ca.gov/waterrights/water_issues/programs/diversion_use/)

coordinates(ag_dat) = c('Longitude','Latitude')#set coordinate system
proj4string(ag_dat)<- "+proj=longlat +ellps=WGS84"
ag_data<-remove.duplicates(ag_dat, remove.second=TRUE) #prioritizes first instance (refer to ordering of joins)
ag_data$use<-"Agriculture"

coordinates(dom_dat) = c('Longitude','Latitude')#set coordinate system
proj4string(dom_dat)<- "+proj=longlat +ellps=WGS84"
dom_data<-remove.duplicates(dom_dat, remove.second=TRUE)
dom_data$use<- "Domestic"

coordinates(ind_dat) = c('Longitude','Latitude')#set coordinate system
proj4string(ind_dat)<- "+proj=longlat +ellps=WGS84"
ind_data<-remove.duplicates(ind_dat, remove.second=TRUE)
ind_data$use <- "Industrial"

coordinates(wild_dat) = c('Longitude','Latitude')#set coordinate system
proj4string(wild_dat)<- "+proj=longlat +ellps=WGS84"
wild_data<-remove.duplicates(wild_dat, remove.second=TRUE)
wild_data$use <- "Fish & Wildlife"

coordinates(other_dat) = c('Longitude','Latitude')#set coordinate system
proj4string(other_dat)<- "+proj=longlat +ellps=WGS84"
other_data<-remove.duplicates(other_dat, remove.second=TRUE)
other_data$use <- "Recreation & Other"

#################################################
###MERGE INDIVIDUAL USE DATA & EXTRACT YEARS#####
#################################################

full.pt<-rbind(ag_data,dom_data,ind_data,wild_data,other_data) #full spatial dataset

# keep only  water rights which had a status change before or during the year of interest
# this keeps out water rights that were activated during a later year 
full.pt.2006 <-full.pt[ full.pt$'Status Date' < as.Date("2007-01-01"),] 
full.pt.2007 <- full.pt[full.pt$'Status Date' < as.Date("2008-01-01"),] 
full.pt.2008 <- full.pt[full.pt$'Status Date'< as.Date("2009-01-01"),] 
full.pt.2009 <- full.pt[full.pt$'Status Date' < as.Date("2010-01-01"),] 
full.pt.2010 <- full.pt[full.pt$'Status Date' < as.Date("2011-01-01"),]
full.pt.2011 <- full.pt[full.pt$'Status Date' < as.Date("2012-01-01"),]
full.pt.2012 <- full.pt[full.pt$'Status Date' < as.Date("2013-01-01"),]
full.pt.2013 <- full.pt[full.pt$'Status Date'< as.Date("2014-01-01"),]
full.pt.2014<- full.pt[full.pt$'Status Date' < as.Date("2015-01-01"),] 
full.pt.2015 <- full.pt[full.pt$'Status Date' < as.Date("2016-01-01"),] 

# remove all water rights which were shut down before the year of interest
# this should leave only water rights which were active for all or some part of the year of interest
notactive <- c("Cancelled", "Closed","Inactive","Rejected","Revoked") 
full.pt.2006<-full.pt.2006[!(full.pt.2006$'Status Date' < as.Date("2006-01-01") & full.pt.2006$POD.Status %in% notactive),] 
full.pt.2007<-full.pt.2007[!(full.pt.2007$'Status Date' < as.Date("2007-01-01") & full.pt.2007$POD.Status %in% notactive),] 
full.pt.2008<-full.pt.2008[!(full.pt.2008$'Status Date' < as.Date("2008-01-01") & full.pt.2008$POD.Status %in% notactive),] 
full.pt.2009<-full.pt.2009[!(full.pt.2009$'Status Date' < as.Date("2009-01-01") & full.pt.2009$POD.Status %in% notactive),] 
full.pt.2010<-full.pt.2010[!(full.pt.2010$'Status Date' < as.Date("2010-01-01") & full.pt.2010$POD.Status %in% notactive),] 
full.pt.2011<-full.pt.2011[!(full.pt.2011$'Status Date' < as.Date("2011-01-01") & full.pt.2011$POD.Status %in% notactive),] 
full.pt.2012<-full.pt.2012[!(full.pt.2012$'Status Date' < as.Date("2012-01-01") & full.pt.2012$POD.Status %in% notactive),] 
full.pt.2013<-full.pt.2013[!(full.pt.2013$'Status Date' < as.Date("2013-01-01") & full.pt.2013$POD.Status %in% notactive),] 
full.pt.2014<-full.pt.2014[!(full.pt.2014$'Status Date' < as.Date("2014-01-01") & full.pt.2014$POD.Status %in% notactive),] 
full.pt.2015<-full.pt.2015[!(full.pt.2015$'Status Date' < as.Date("2015-01-01") & full.pt.2015$POD.Status %in% notactive),] 

######################################################
###WRITE WATER RIGHT POINT DATA to ESRI shapefiles ###
######################################################

writeOGR(full.pt.2006,".","full.pt.2006", driver="ESRI Shapefile")
writeOGR(full.pt.2007,".","full.pt.2007", driver="ESRI Shapefile")
writeOGR(full.pt.2008,".","full.pt.2008", driver="ESRI Shapefile")
writeOGR(full.pt.2009,".","full.pt.2009", driver="ESRI Shapefile")
writeOGR(full.pt.2010,".","full.pt.2010", driver="ESRI Shapefile")
writeOGR(full.pt.2011,".","full.pt.2011", driver="ESRI Shapefile")
writeOGR(full.pt.2012,".","full.pt.2012", driver="ESRI Shapefile")
writeOGR(full.pt.2013,".","full.pt.2013", driver="ESRI Shapefile")
writeOGR(full.pt.2014,".","full.pt.2014", driver="ESRI Shapefile")
writeOGR(full.pt.2015,".","full.pt.2015", driver="ESRI Shapefile")


