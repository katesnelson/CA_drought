library(nlme)
library(dplyr)
library (lme4)
library(rgdal)
library(INLA)
inla.upgrade(testing = T)
library(spdep)
library(maptools)

#read in data
dsf <- tbl_df(read.csv("C:\\Users\\Emily Burchfield\\delta_lcopbox\\Vanderbilt\\Kate_Emily\\Data\\pt_data_kn.txt", stringsAsFactors = FALSE))
ds = dsf[dsf$farm_flag==1,]

#change datatype 
ds = transform(ds, 
               tvp_09 = as.numeric(tvp_09),
               tvp_14 = as.numeric(tvp_14),
               delta_tvp = as.numeric(delta_tvp),
               AreaAcres = as.numeric(AreaAcres))

#change values of 999.99 to missing value (Perc_Rip and Perc_Pre1914)
ds[ds==999.99] <- NA

#farming pixels only
data <- ds[ds$farm_flag == 1,]

Rip <- ds$COUNT_Rip/ds$fmmp_area
P1914 <- ds$COUNT_Pre1914/ds$fmmp_area

#read in shapefile
#https://rpubs.com/corey_sparks/132760
coord <- readOGR("C:\\Users\\Emily Burchfield\\delta_lcopbox\\Vanderbilt\\Kate_Emily\\Data", "HUC12")













#playing with L2
coord <- readOGR("C:\\Users\\Emily Burchfield\\delta_lcopbox\\Vanderbilt\\Kate_Emily\\Data", "HUC_centroid")
x <- coordinates(coord)[,-2] #lon
y <- coordinates(coord)[,-1] #lat
sp1 <- corSpatial(1, form = ~x + y, type = "gaussian")
scor <- Initialize(sp1, as(coord, "data.frame") [,c("lon", "lat")], nugget = F)

M4 <- lmer(delta_tvp ~ delta_lc + diverse + WR_density + P1914 + GW_dnsty15 + 
             delta_lc*diverse + (1 + delta_lc|dsf), data = ds, REML = "False", correlation = scor) 
summary(M4)




#https://HUC12s.google.com/forum/#!topic/r-inla-discussion-HUC12/DYcAjPxlAz8

# random intercept
m1 = lmer(delta_lc ~ GW_dnsty15 + (1|HUC12), data = dsf)	
form1 = delta_lc ~ GW_dnsty15 + f(HUC12, model = 'iid')

#organize data
n.HUC12 = length(dsf$HUC12) #make sure this is what he was getting at here
dsf$i.int = dsf$HUC12
dsf$j.int = dsf$HUC12 + n.HUC12
dsf$k.int = dsf$j.int + n.HUC12

# uncorrelated random intercept and random slope
m2 = lmer(delta_lc ~ GW_dnsty15 + (1|HUC12) + (0 + GW_dnsty15|HUC12), data = dsf)
form2 = delta_lc ~ GW_dnsty15 + f(i.int, model = 'iid') + f(j.int, GW_dnsty15, model = 'iid')

# correlated random intercept and random slope
m3 = lmer(delta_lc ~ GW_dnsty15 + (GW_dnsty15 |HUC12), data = dsf)
form3 = delta_lc ~ GW_dnsty15 + f(i.int, model = 'iid2d', n = 2*n.HUC12) + f(j.int, GW_dnsty15, copy = 'i.int')

# uncorrelated random intercept and random slope of GW_dnsty15^2
m4 = lmer(delta_lc ~ GW_dnsty15 + I(GW_dnsty15^2) + (1|HUC12) + (0 + I(GW_dnsty15^2)|HUC12), data = dsf)
form4 = delta_lc ~ GW_dnsty15 + I(GW_dnsty15^2) + f(i.int, model = 'iid') + f(j.int, I(GW_dnsty15^2), model = 'iid')


#random intercept and slopes, with correlated slopes
m5 = lmer(delta_lc ~ GW_dnsty15 + I(GW_dnsty15^2) + (1|HUC12) + (0 + GW_dnsty15 + I(GW_dnsty15^2)|HUC12), data = dsf)
form = delta_lc ~ GW_dnsty15 + I(GW_dnsty15^2) + f(HUC12, model = 'iid') + f(i.int, GW_dnsty15, model = 'iid2d', n = 2*n.HUC12) + f(j.int, I(GW_dnsty15^2), copy = 'i.int')


# correlated random intercept and slopes
m6 = lmer(delta_lc ~ GW_dnsty15 + I(GW_dnsty15^2) + (GW_dnsty15 + I(GW_dnsty15^2)|HUC12), data = dsf)
form5 = delta_lc ~ GW_dnsty15 + I(GW_dnsty15^2) + f(i.int, model = 'iid3d', n = 3*n.HUC12) + f(j.int, GW_dnsty15, copy = 'i.int') + f(k.int, I(GW_dnsty15^2), copy = 'i.int')





