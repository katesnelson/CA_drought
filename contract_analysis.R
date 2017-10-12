library(rgdal)
library(pdftools)

#This script searches lists of agriculutral water contractors in California for matches with owners of water rights in California. Looser matches between water right
#owner names and listed water contractors were discovered by manually searching the list of water rights owners. All indentified matches in the water rights dataset 
#were mapped to examine the spatial distribution of water rights PODs associated with contract water. Identified matches by water project source are provided at the 
#end of this script.


# download the file
download.file("https://www.usbr.gov/mp/cvp-water/docs/ag-contractors-website.pdf", ".\\data\\contract.pdf", mode="wb")
txt <- pdf_text(".\\data\\contract.pdf")

# get shapefile data
pt <- readOGR(dsn = "C:\\Users\\Emily\\Dropbox\\Vanderbilt\\Kate_Emily\\Data\\full.pt.2015.shp")
# subset agriculture PODs
pt_ag <- pt[pt$use == "Agriculture",]
owner <- as.character(unique(pt_ag$PrmryOw))

#clean up txt file
txt_split <- unlist(strsplit(txt, "\r\n"))
txt_split <- unlist(strsplit(txt_split, "   "))
txt_split <- txt_split[txt_split != ""]
txt_sub <- txt_split[which(txt_split != txt_split[2])] # drops No.
txt_sub <- txt_sub[-grep("Contract*", txt_sub)]
txt_sub <- txt_sub[-grep("Friday*", txt_sub)]
txt_sub <- txt_sub[-grep("Page*", txt_sub)]
txt_sub <- txt_sub[-grep("^\\s+[0-9]+$", txt_sub)]
contracts <- trimws(txt_sub, which="both")


#https://www.r-bloggers.com/fuzzy-string-matching-a-survival-skill-to-tackle-unstructured-information/

#compute distance
dist.name <- adist(owner, contracts, ignore.case=T)

#take pairs with min distance
min.name <- apply(dist.name, 1, min)

match.s1.s2<-NULL  
for(i in 1:nrow(dist.name))
{
  s2.i<-match(min.name[i],dist.name[i,])
  s1.i<-i
  match.s1.s2<-rbind(data.frame(s2.i=s2.i,s1.i=s1.i,s2name=contracts[s2.i], s1name=owner[s1.i], adist=min.name[i]),match.s1.s2)
}
# and we then can have a look at the results
match <- match.s1.s2
colnames(match) <- c("s2", "s1", "contract_name", "owner_name", "adist")

# seems like adist < 3 shows actual matches
match <- match[match$adist < 3,]

#drop NA
match <- match[-which(is.na(match)),]

# 215 contracters listed on the PDF.  Here, we find a match of only 30 contracters
match$owner_name

# add the state water project contractors:
# http://www.water.ca.gov/swpao/wsc.cfm

# 32 listed on website

# didn't find: Antelope, East Kern, Castaic, Crestline-Lake Arrowhead,Devil's Den,
# dudley, empire west, hacienda, kings county, mojave, napa, oak flat, san gabriel,
# san gornio, san luis, tulare, west covina, yuba

swc <- c(owner[grep("ALAMEDA COUNTY WATER DISTRICT*", owner, ignore.case=T)],
         owner[grep("ALAMEDA COUNTY F C & W C D (ZONE 7)*", owner, ignore.case=T)],
         owner[grep("Butte Valley Irrigation*", owner, ignore.case=T)], 
         owner[grep("Coachella*", owner, ignore.case=T)], 
         owner[grep("kern county*", owner, ignore.case=T)],
         owner[grep("littlerock*", owner, ignore.case=T)],
         owner[grep("metropol*", owner, ignore.case=T)],
         owner[grep("palmdale*", owner, ignore.case=T)], 
         owner[grep("plumas county*", owner, ignore.case=T)],
         owner[grep("san bernardino valley*", owner, ignore.case=T)],
         owner[grep("santa barbara county*", owner, ignore.case=T)],
         owner[grep("santa clara valley water*", owner, ignore.case=T)],
         owner[grep("solano*", owner, ignore.case=T)],
         owner[grep("ventura county watershed*", owner, ignore.case=T)])

# irrigation district == water district? (palmdale, plumas county)

################################################################################################################################################
#MAP identifed matches

d<-readRDS('cv.dat06072017.rds')
full<-shapefile('C:/Users/tuan/Dropbox/Kate_Emily/Data/full.pt.2015.shp')
cv <- shapefile('cv_huc.shp')

contractnames<-c("STATE OF CALIFORNIA - DEPARTMENT OF WATER RESOU", "CALIF DEPT OF WATER RESOURCES", 
                 "STATE WATER RESOURCES CONTROL BOARD","STOCKTON EAST WATER DISTRICT", "KINGS RIVER WATER ASSOCIATION", 
                 "???KERN County Water Agency",  "KERN DELTA WATER DISTRICT", "Buena Vista Water Storage District",
                 "CAWELO WATER DISTRICT", "ROSEDALE- RIO BRAVO WATER STORAGE DISTRICT",  "SEMITROPIC WATER STORAGE DISTRICT",
                 "A L ANDERSON", "ANDERSON-COTTONWOOD IRRIGATION DISTRICT", "Arnold, Mike, Mark Andreotti, and Assoc L",
                 "Banta-Carbona Irrigation District", "CHRISTOPHER BARDIS", "DIANNE E BUTLER", "Byrd, Ann & Osborne, Jane",
                 "BYRON-BETHANY IRRIGATION DISTRICT", "Cachil Dehe Band of Wintun Indians", "CARTER MUTUAL WATER COMPANY",
                 "CHOWCHILLA WATER DISTRICT", "Mary Coelho Trust", "CONAWAY PRESERVATION GROUP LL", "DRISCOLL STRAWBERRY ASSOCIATES",
                 "EASTSIDE MUTUAL WATER COMPANY",  "WALLACE L EDSON", "ALLEN A EHRKE", "FEATHER WATER DISTRICT", "FEDORA FARMS INC",
                 "FORESTHILL PUBLIC UTILITY DISTRICT", "GLENN-COLUSA IRRIGATION DISTRICT", "GRAVELLY FORD WATER DISTRICT",               
                 "MILDRED HEIDRICK", "Heidrick & Heidrick Properties LP", "HENLE FAMILY LIMITED PARTNERSHIP",
                 "HERSHEY LAND COMPANY ROW CROP LLC", "HOWALD FARMS", "CAROL KARY", "KNAGGS FARMING COMPANY LP",
                 "WILLIAM P LOCKETT", "LOMO COLD STORAGE A CAL GEN PARTNERSHIP", "Lonon, Michael E. - MICHAEL E LONON",
                 "LOWER TULE RIVER IRRIGATION DISTRICT", "MADERA IRRIGATION DISTRICT", "MAXWELL IRRIGATION DISTRICT",                               
                 "MERIDIAN FARMS WATER CO", "JOSEPH A MOREHEAD", "MUNSON FAMILY TRUST", "NATOMAS CENTRAL MUTUAL WATER CO",
                 "THOMAS L NELSON & HAZEL M NELSON TRUST", "Frank J. O'Brien Family Trust", "ODYSSEUS FARMS",                              
                 "Oji Brothers Farm, Inc.", "ORANGE COVE IRRIGATION DISTRICT", "PAJARO VALLEY WATER MANAGEMENT AGENCY",
                 "SANTA CLARA VALLEY WATER DISTRICT", "PATTERSON IRRIGATION DISTRICT", "PELGER MUTUAL WATER COMPANY",                
                 "Pixley Irrigation District", "PLACER COUNTY WATER AGENCY", "PRINCETON-CODORA-GLENN IRRIGATION DISTRICT",
                 "PROVIDENT IRRIGATION DISTRICT", "Reclamation District No. 1004", "RECLAMATION DISTRICT #1004",
                 "RECLAMATION DISTRICT #108", "HENRY D RICHTER JR", "RIVER GARDEN FARMS COMPANY", "SAN BENITO COUNTY WATER DISTRICT",          
                 "SANTA CLARA VALLEY WATER DISTRICT", "CHARLES W SEAVER", "SIDDIQUI FAMILY PARTNERSHIP", "SIOUX CREEK PROPERTY LL",         
                 "SUTTER MUTUAL WATER COMPANY", "STONY CREEK WATER DISTRICT", "STOCKTON EAST WATER DISTRICT", "STEVE TARKE",
                 "WEST SIDE IRRIGATION DISTRICT", "Tisdale Irrigation & Drainage CO", "Thomes Creek Water Users Association", 
                 "Sycamore Family Trust", "TERRA BELLA IRRIGATION DISTRICT", "TULARE IRRIGATION DISTRICT", "KENNETH L WALLACE",
                 "WEST STANISLAUS IRRIGATION DISTRICT", "WINDSWEPT LAND AND CATTLE COMPANY", "ALAMEDA COUNTY WATER DISTRICT",                         
                 "ALAMEDA COUNTY F C & W C D (ZONE 7)", "BUTTE VALLEY IRRIGATION DISTRICT", "COACHELLA VALLEY WATER DISTRICT",                       
                 "KERN COUNTY WATER AGENCY", "LITTLEROCK CREEK IRRIGATION DISTRICT", "The Metropolitan Water District of Southern California",
                 "PALMDALE IRRIGATION DISTRICT", "PLUMAS COUNTY", "SAN BERNARDINO VALLEY W C D", "SANTA BARBARA COUNTY F C and W C DISTRICT",             
                 "SANTA CLARA VALLEY WATER DISTRICT", "SOLANO COUNTY WATER AGENCY", "VENTURA COUNTY WATERSHED PROTECTION DISTRICT",
                 "DESERT WATER AGENCY", "COUNTY OF SAN LUIS OBISPO")

sub<-full[full@data$PrmryOw %in% contractnames, ]
sub <- spTransform(sub, CRS(proj4string(cv)))
subcv<-sub[cv, ]

#dsub<-d@data[d@data$flflag==1,"GRID_CODE" ]
#cv2<-cv[d]

plot(cv)
points(subcv, cex=0.8, bg="red", col="red", pch=19)

library(ggmap)
library(RColorBrewer) 

colors <- brewer.pal(9, "RdPu")

mapImage <- get_map(location = c(lon = -120, lat = 37.5),
                    color = "color",
                    source = "google",
                    maptype = "terrain",
                    zoom = 6)

cv.map <- spTransform(cv, '+init=epsg:4326 +proj=longlat
                      +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

subcv <- spTransform(subcv, '+init=epsg:4326 +proj=longlat
                      +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

plot(mapImage)
plot(cv.map, add=TRUE)
points(subcv, cex=0.8, bg="red", col="red", pch=19)

##########################################################################################################################################################
#All matches, identified using text analysis and manual searches allowing for multiple spellings/versions of names and association of multiple 
# associations (subcontractors) within large water districts with a single contractor listing

# Central Valley Project (CVP) Water Contractors with USBR
# 
# 20 	California, State of 
# ???  STATE OF CALIFORNIA - DEPARTMENT OF WATER RESOU
# ??? CALIF DEPT OF WATER RESOURCES
# ??? STATE WATER RESOURCES CONTROL BOARD
# 
# 22 	Central San Joaquin Water Conservation District
# ??? STOCKTON EAST WATER DISTRICT (http://csjwcd.com/about-csjwcd/history/)
# 
# 60 	Fresno Irrigation District 
# ??? KINGS RIVER WATER ASSOCIATION (http://www.krcd.org/_pdf/Kings_River_Handbook_2009.pdf)
# 
# 96 	Kern-Tulare Water District (also called the Water Association of Kern County, http://www.wakc.com/whos-who/kern-tulare-water-district/)
# ??? KERN County Water Agency (http://www.wakc.com/wp-content/uploads/2016/01/SWP-Contracts-in-Kern-County.pdf)
# ??? KERN DELTA WATER DISTRICT
# ??? Buena Vista Water Storage District
# ??? CAWELO WATER DISTRICT
# ??? ROSEDALE- RIO BRAVO WATER STORAGE DISTRICT
# ??? SEMITROPIC WATER STORAGE DISTRICT
# 
# 
# 4 	A L ANDERSON
# 6 	ANDERSON-COTTONWOOD IRRIGATION DISTRICT   
# 7 	Arnold, Mike, Mark Andreotti, and Assoc L
# 10 	Banta-Carbona Irrigation District          
# 11 	CHRISTOPHER BARDIS
# 15 	DIANNE E BUTLER
# 17 	Byrd, Ann & Osborne, Jane
# 18 	BYRON-BETHANY IRRIGATION DISTRICT  
# 19 	Cachil Dehe Band of Wintun Indians         
# 21 	CARTER MUTUAL WATER COMPANY
# 25 	CHOWCHILLA WATER DISTRICT                 
# 28 	Mary Coelho Trust
# 31 	CONAWAY PRESERVATION GROUP LL
# 42 	DRISCOLL STRAWBERRY ASSOCIATES
# 49  	EASTSIDE MUTUAL WATER COMPANY      
# 50 	WALLACE L EDSON
# 52 	ALLEN A EHRKE
# 55 	FEATHER WATER DISTRICT                     
# 56 	FEDORA FARMS INC
# 57  	FORESTHILL PUBLIC UTILITY DISTRICT         
# 70   	GLENN-COLUSA IRRIGATION DISTRICT           
# 73    GRAVELLY FORD WATER DISTRICT               
# 78	MILDRED HEIDRICK 
# 77 	Heidrick & Heidrick Properties LP
# 79 	HENLE FAMILY LIMITED PARTNERSHIP
# 80 	HERSHEY LAND COMPANY ROW CROP LLC
# 85 	HOWALD FARMS
# 95 	CAROL KARY
# 101 	KNAGGS FARMING COMPANY LP
# 110 	WILLIAM P LOCKETT
# 112 	LOMO COLD STORAGE A CAL GEN PARTNERSHIP
# 113 	Lonon, Michael E. - MICHAEL E LONON
# 114   LOWER TULE RIVER IRRIGATION DISTRICT      
# 116 	MADERA IRRIGATION DISTRICT
# 117 	MAXWELL IRRIGATION DISTRICT                               
# 122 	MERIDIAN FARMS WATER CO
# 124 	JOSEPH A MOREHEAD
# 126	MUNSON FAMILY TRUST
# 129	NATOMAS CENTRAL MUTUAL WATER CO
# 130	THOMAS L NELSON & HAZEL M NELSON TRUST
# 132	Frank J. O'Brien Family Trust
# 133	ODYSSEUS FARMS                             
# 135	Oji Brothers Farm, Inc.  
# 137	ORANGE COVE IRRIGATION DISTRICT
# 143	PAJARO VALLEY WATER MANAGEMENT AGENCY
# 143	SANTA CLARA VALLEY WATER DISTRICT
# 145   PATTERSON IRRIGATION DISTRICT             
# 146	PELGER MUTUAL WATER COMPANY                
# 149   Pixley Irrigation District                
# 150	PLACER COUNTY WATER AGENCY                
# 153	PRINCETON-CODORA-GLENN IRRIGATION DISTRICT
# 155	PROVIDENT IRRIGATION DISTRICT             
# 160	Reclamation District No. 1004 
# 160	RECLAMATION DISTRICT #1004
# 161	RECLAMATION DISTRICT #108
# 166	HENRY D RICHTER JR            
# 167	RIVER GARDEN FARMS COMPANY                
# 173	SAN BENITO COUNTY WATER DISTRICT          
# 175	SANTA CLARA VALLEY WATER DISTRICT 
# 177	CHARLES W SEAVER
# 181	SIDDIQUI FAMILY PARTNERSHIP
# 182	SIOUX CREEK PROPERTY LL         
# 189	SUTTER MUTUAL WATER COMPANY               
# 188	STONY CREEK WATER DISTRICT                 
# 186	STOCKTON EAST WATER DISTRICT
# 191	STEVE TARKE 
# 195	WEST SIDE IRRIGATION DISTRICT
# 197 	Tisdale Irrigation & Drainage CO
# 196 	Thomes Creek Water Users Association
# 190	Sycamore Family Trust
# 194 	TERRA BELLA IRRIGATION DISTRICT
# 201	TULARE IRRIGATION DISTRICT 
# 205	KENNETH L WALLACE
# 206	WEST STANISLAUS IRRIGATION DISTRICT   
# 211	WINDSWEPT LAND AND CATTLE COMPANY                     
# 
# 
# 
# State Water Project-Water Supply Contracts and Amendments
# 
# [1] "ALAMEDA COUNTY WATER DISTRICT"                         
# [2] "ALAMEDA COUNTY F C & W C D (ZONE 7)"                   
# [3] "BUTTE VALLEY IRRIGATION DISTRICT"                      
# [4] "COACHELLA VALLEY WATER DISTRICT"                       
# [5] "KERN COUNTY WATER AGENCY"                              
# [6] "LITTLEROCK CREEK IRRIGATION DISTRICT"                  
# [7] "The Metropolitan Water District of Southern California"
# [8] "PALMDALE IRRIGATION DISTRICT"                          
# [9] "PLUMAS COUNTY"                                         
# [10] "SAN BERNARDINO VALLEY W C D"                           
# [11] "SANTA BARBARA COUNTY F C and W C DISTRICT"             
# [12] "SANTA CLARA VALLEY WATER DISTRICT"                     
# [13] "SOLANO COUNTY WATER AGENCY"                            
# [14] "VENTURA COUNTY WATERSHED PROTECTION DISTRICT"
# DESERT WATER AGENCY
# COUNTY OF SAN LUIS OBISPO

