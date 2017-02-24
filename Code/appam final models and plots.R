library(memisc)
library(dplyr)
library (gstat)
library(sp)
library(nlme)
library(rgdal)
library(INLA)
library(spdep)
library(raster)


setwd('C:/Users/nelsonks/Desktop/drought')
setwd('E:/drought')
#setwd('/home/kate/R/x86_64-pc-linux-gnu-library/3.3/INLA/')
##############################################
#Read in datasets to use for inla models
##############################################
#full cleaned pixel dataset of central valley as spatial points data frame
#d <- readRDS('/data/emily/WF/kate/final_data/cv.dat.clean.rds')

###  full dataset with lat lon columns and updated wr info --> check for cleanliness 
#d<- readRDS('/home/kate/CA/data/cv.dat2.rds') 
d<- readRDS('cv.dat2.rds') 
#plot(d)


#full central valley watershed as spatial polygons data frame
#cv <- shapefile('/data/emily/WF/kate/shp/cv_huc.shp')
cv <- shapefile('cv_huc.shp')
cv <- spTransform(cv, CRS(proj4string(d))) #set to same proj
proj4string(cv)==proj4string(d) #double check proj

#plot together to check validity of spatial overlay
#plot(d, cex=0.05, col="blue")
#lines(cv)

#remove uneeded info from cv
cv<-cv[ ,8]



##############################################
#Select farmland subset for analysis
##############################################
#select pixel record dataframe based on farmland flag

ds <- d@data[d@data$flflag ==1, ]

#select subset of cv for creating adjancency matrix
list<-ds$huc12_id
cv$HUC12<-as.numeric(cv$HUC12)
cvs<-cv[cv@data$HUC12 %in% list, ]
#plot(cvs)

################################################
#Create adjacency matrix
############################################

shp<-cvs

#order ascending by HUC
shp <- shp[order(shp$HUC12),]

####build group index that corresponds to adjacency matrix#####
HUC<- seq(1,length(shp@data[ ,1]),1)
shp@data<- cbind(shp@data, HUC)
#plot(shp)

#join new simple index to full ds####
ds<-left_join(ds,shp@data, by=c("huc12_id"= "HUC12"))

#https://groups.google.com/forum/#!topic/r-inla-discussion-group/iQcOqUHwxjM
# ds <- ds[order(ds$year,ds$HUC),] #order by time then huc
#try ordering by huc then time like they do in the Ohio example chapter 7 of r-inla book --> works for bym model
ds <- ds[order(ds$HUC,ds$year),]

#create neighbors list from polygon object
temp <- poly2nb(shp, queen=FALSE) #neighbors must share more than one point at boundary

#convert to a sparse matrix to reduce memory
H.adj <- nb2mat(temp, style ="B") #convert neighbor list to binary coded neighbor weights matrix
H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse sytle matrix

#plot spatial object
image(inla.graph2matrix(H.adj), xlab="", ylab="")

# #write out list as an INLA object to wd
# nb2INLA("/home/kate/CA/data/H.adj", H.adj)

#############################################
#Create SPDE object
##############################################
coords <- ds[ ,c(4,3)] #create a list of coordinates from lat,lon in ds
coords<-unique(coords) # use unique coords to create mesh
mesh<-inla.mesh.2d(loc=coords,  max.edge=200000, cutoff = 5000 ) #create a mesh for estimation, trouble runing when specify an offset
plot(mesh) #check projection stuff
spde <-inla.spde2.matern(mesh=mesh, alpha=2) #pg 206, this mesh is for accounting for continuous spatial processes (landscape,etc..), and not for interpolating predictor values
# creating the observation/projection matrix (sparse weight matrix pg 205)

coords.alltime<-as.matrix(coords[ds$POINTID, c("lon", "lat")]) #list of coordinates for all times
#y<-ds$year + 1
#n.year<-length(unique(y))
A.est <-inla.spde.make.A(mesh=mesh,loc=coords.alltime)#, group=y, n.group=n.year) #object that reconciles mesh and observations

s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)

########################################
#Prep data for analysis
#######################################

y=ds$px_tvp #dependent variable


ds$year <- ds$year - 7 #center time at 2007 = 0
ds <- mutate(ds, year1=year, year2=year, year3=year, year4 = year) #add duplicate year grouping index
ds<-mutate(ds,POINTID1=POINTID)#add duplicate pixel grouping index
ds <- mutate(ds, HUC1=HUC,HUC2=HUC, HUC3=HUC, HUC4=HUC, HUC5=HUC, HUC6=HUC, HUC7=HUC, HUC8=HUC,HUC9=HUC, 
             HUC10=HUC, HUC11=HUC, HUC12=HUC, HUC13=HUC, HUC14=HUC, HUC15=HUC, HUC16=HUC, HUC17=HUC) #add duplicates of HUC index for random slope grouping
ds<- mutate(ds, wr_dens1=wr_dens, rip_perc1=rip_perc, approp_perc1=approp_perc, pre_perc1=pre_perc)
s.index1 <- inla.spde.make.index(name="spatial.field1", n.spde=spde$n.spde)

#standardize independent vars to ease interpretation
ds$gw<- scale(ds$gw)
#ds$spi<-scale(ds$spi) #already on a z-score scale
ds$fl_cd<-scale(ds$fl_cd)
ds$wr_dens<-scale(ds$wr_dens)
#ds<-mutate(ds,ratio=(rip_cnt+1)/(approp_cnt +1))


##################################################
#Models for appam
##################################################


##Null 3 level MLM ## 3 level model for instances (time) nested in pixels nested in hucs, use to get 3 level ICC --1.5 hr run time
null_3lvl2  <-y ~ 1 + f(HUC, model="iid") + f(POINTID, model = "iid")
output.null.3lvl2 <- inla (null_3lvl2, family = "gaussian", data=ds, verbose=TRUE)
saveRDS(output.null.3lvl2, file='/home/kate/CA/data/output.null.3lvl3')
n33<-readRDS('/home/kate/CA/data/output.null.3lvl3')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm3c<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens + rip_perc + approp_perc + pre_perc +
  year*rip_perc +year*approp_perc + year*pre_perc +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm3c <- inla (stlm3c,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm3c, file='/home/kate/CA/data/output.stlm3c')
stlm3c<-readRDS('output.stlm3c')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm3d<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens + rip_perc + approp_perc + pre_perc +
  spi*rip_perc +spi*approp_perc + spi*pre_perc +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm3d <- inla (stlm3d,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm3d, file='/home/kate/CA/data/output.stlm3d')
stlm3d<-readRDS('/home/kate/CA/data/output.stlm3d')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC and POINTID, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm3e<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens + rip_perc + approp_perc + pre_perc +
  year*rip_perc +year*approp_perc + year*pre_perc +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1") + f(POINTID, model="iid")
output.stlm3e <- inla (stlm3e,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm3e, file='/home/kate/CA/data/output.stlm3e')
stlm3e<-readRDS('/home/kate/CA/data/output.stlm3e')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm4a<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens +  approp_perc + pre_perc +
  year*approp_perc + year*pre_perc +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm4a <- inla (stlm4a,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm4a, file='output.stlm4a')
stlm4a<-readRDS('output.stlm4a')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm4b<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens +  rip_perc + pre_perc +
  year*rip_perc + year*pre_perc +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm4b <- inla (stlm4b,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm4b, file='output.stlm4b')
stlm4b<-readRDS('output.stlm4b')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm4c<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens +  approp_perc + rip_perc +
  year*approp_perc + year*rip_perc +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm4c <- inla (stlm4c,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm4c, file='output.stlm4c')
stlm4c<-readRDS('output.stlm4c')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm4d<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens +  approp_perc  +
  year*approp_perc  +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm4d <- inla (stlm4d,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm4d, file='output.stlm4d')
stlm4d<-readRDS('output.stlm4d')

#  model with spatial random effects and year fixed effect and space-time random effects(global time effect, random intercepts by HUC, predictor-time interactions, and differential time trend by HUC, and ar1 time error term)
stlm4e<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens +  pre_perc  +
  year*pre_perc  +
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
output.stlm4e <- inla (stlm4e,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.stlm4e, file='output.stlm4e')
stlm4e<-readRDS('output.stlm4e')
###################################
#Variance Extracts for null model and ICC
#########################

prec<-n$marginals.hyperpar$`Precision for the Gaussian observations`
eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-n$marginals.hyperpar$`Precision for HUC`
u00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2)
prec3<-n$marginals.hyperpar$`Precision for POINTID`
u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)

ICC1<-(eijk/(eijk+u0jk+u00k))
ICC2<-(u0jk/(u0jk+u00k+eijk))
ICC3<-(u00k/(u00k+u0jk+eijk))


###############################
#Maps & Plots
#############################

##FIGURE 1
#central valley spatail extents map
library(ggmap)
library(RColorBrewer) 

colors <- brewer.pal(9, "RdPu")

mapImage <- get_map(location = c(lon = -120, lat = 37.5),
                    color = "color",
                    source = "google",
                    maptype = "terrain",
                    zoom = 6)

cv.map <- spTransform(cvs, '+init=epsg:4326 +proj=longlat
                      +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
cv.map<-fortify(cv.map)

ggmap(mapImage) +
  geom_polygon(aes(x = long,
                   y = lat,
                   group = group),
               data = cv.map,
               color = colors[9],
               fill = colors[6],
               alpha = 0.3) +
  labs(x = "Longitude",
       y = "Latitude")

#################################################
#KRIG Plots
##################################################
pred3<-readRDS('krigged_data.rds')

stplot(pred3[,c(85,91,97,103,109,121,133,145,157,169)], main="Interpolated Groundwater Elevations (feet above mean sea level)", cex.main=1)

kd<-as(pred3, "data.frame") #transform the stfdf to a data frame

pred_time_trend <- kd%>% group_by(time) %>% summarize(mean(var1.pred)) 
plot(x=pred_time_trend$time, y=pred_time_trend$`mean(var1.pred)`) #plot mean of predicted elev across space over time

pred_space_trend <- kd%>% group_by(sp.ID) %>% summarize(mean=mean(var1.pred)) 
ksp<- left_join(kd, pred_space_trend, by = "sp.ID")
coordinates(ksp)<- ~x1-x2
ksp@proj4string <- CRS('+init=epsg:3310')
gridded(ksp)=TRUE
spplot(ksp[ ,6], main="Average Predicted Groundwater Elevation \n(ft above mean sea level) for 2007-2014", cex.main=1) #map mean of predicted elev across time over space

########################################

sub<-d
sub<-sub[sub@data$flflag==1,]
sub@data<-sub@data[sub@data$year==14,] #pull out one year at a time
spplot(sub[ ,"AreaSqKm"], cex=0.5)
spplot(sub[ ,"wr_cnt"], cex=0.5)
spplot(sub[ ,"wr_dens"], cex=0.5)
spplot(sub[ ,"pre_cnt"], cex=0.5)
spplot(sub[ , "ag_cnt"], cex=0.5)
spplot(sub[ , "adjud_cnt"], cex=0.5)
spplot(sub[ , "gw"], cex=0.5) #looks good
spplot(sub[ ,"diverse"])
spplot(sub[ ,"fl_cd"])
spplot(sub[ ,"adjud_perc"])
spplot(sub[ ,"pre_perc"])
spplot(sub[ ,"px_tvp"])
spplot(sub[sub@data$px_tvp==0, ])


##FIGURE 2a
hist(ds$px_tvp, main="Histogram of All TVP Values", xlab="TVP")

##FIGURE 2b
#tvp mean trend through time
t<-ds %>% group_by(year) %>%summarise(mean=mean(px_tvp))
t2<-ds %>% group_by(year) %>%summarise(sd=sd(px_tvp))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim = c(0,1), main="Mean TVP through Time", xlab="Years since 2007", ylab="Mean TVP for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

#tvp mean trend through space
shpdup<-shp 
s<-ds %>% group_by(HUC) %>%summarise(mean(px_tvp))
s2<-ds %>% group_by(HUC) %>%summarise(sd(px_tvp))
shp@data<-full_join(shp@data, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","tvp")
spplot(shp[ ,3])

#tvp mean trends through space time
st<-ds[ds$year==0, ]
st1<- st%>%group_by(HUC) %>% summarise(mean(px_tvp))
shp@data<-full_join(shp@data, st1, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","tvp", "tvp0")
spplot(shp[ ,4])

st<-ds[ds$year==7, ]
st2<- st%>%group_by(HUC) %>% summarise(mean(px_tvp))
shp@data<-full_join(shp@data, st2, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","tvp", "TVP_2007", "TVP_2014")
spplot(shp[ ,5])

shp@data<-mutate(shp@data, Delta_TVP=(tvp7 - tvp0))
spplot(shp[ ,6])

##FIGURE 3
spplot(shp,c("TVP_2007","TVP_2014"), main="Mean TVP for each Watershed")

#FIGURE 4
sub<-d
sub@data<-sub@data[sub@data$year==14,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]#pull the farmland area
z1<- spplot(sub[, "px_tvp"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="TVP in 2014", cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z1$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z1)
sub<-d
sub@data<-sub@data[sub@data$year==7,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]
z2<- spplot(sub[, "px_tvp"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),  main="TVP in 2007", cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z2$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z2)
print(z2, position=c(0,0,0.5,1), more=TRUE)
print(z1, position=c(0.5,0,1,1))

#FIGURE 9 --> hmm why so many negaitive values? (how was gw calculated?) (end of yr - beg of yr)
sub<-d
sub@data<-sub@data[sub@data$year==14,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]#pull the farmland area
z1<- spplot(sub[, "gw"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="TVP in 2014")#, cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z1$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z1)
sub<-d
sub@data<-sub@data[sub@data$year==7,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]
z2<- spplot(sub[, "gw"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),  main="TVP in 2007")#, cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z2$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z2)
print(z2, position=c(0,0,0.5,1), more=TRUE)
print(z1, position=c(0.5,0,1,1))

####################

 #attr from scaling throws off the summarise use non-scaled data
#FIGURE XX
#spi mean trend through time
t<-ds %>% group_by(year) %>% summarise(mean=mean(spi))
t2<-ds %>% group_by(year) %>%summarise(sd=sd(spi))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim = c(-7,4), main="Mean SPI through Time", xlab="Years since 2007", ylab="Mean SPI for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

#spi mean trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(spi))
s2<-ds %>% group_by(HUC) %>%summarise(sd(spi))
shp2<-shp[ ,1:2]
shp2@data<-full_join(shp2@data, s, by="HUC")
colnames(shp2@data)<-c("HUC12", "HUC","spi")
spplot(shp2[ ,3])

#spi mean trends through space time
st<-ds[ds$year==0, ]
st1<- st%>%group_by(HUC) %>% summarise(mean(spi))
shp2@data<-full_join(shp2@data, st1, by="HUC")
colnames(shp2@data)<-c("HUC12", "HUC","spi", "spi2007")
spplot(shp2[ ,4])

st<-ds[ds$year==7, ]
st2<- st%>%group_by(HUC) %>% summarise(mean(spi))
shp2@data<-full_join(shp2@data, st2, by="HUC")
colnames(shp2@data)<-c("HUC12", "HUC","SPI", "SPI_2007", "SPI_2014")

###FIGUREXX
spplot(shp2, c("SPI_2007","SPI_2014"), main ="Mean SPI for each Watershed")


###These are pretty much flat###

# #wr_dens mean trend through time
# wr<-ds %>% group_by(year) %>%summarise(mean((wr_dens)))
# wr2<-ds %>% group_by(year) %>%summarise(sd((wr_dens)))
# plot(wr, ylim=c(0,1))
# 
# #rip mean trend through time
# r<-ds %>% group_by(year) %>%summarise(median(rip_perc)) #problem with na and inf
# r2<-ds %>% group_by(year) %>%summarise(sd(!is.na(rip_perc)))
# points(ds$year, ds$rip_perc, ylim=c(0.6,1),col="red", cex=0.5)
# 
# #approp mean trend through time
# a<-ds %>% group_by(year) %>%summarise(mean(!is.na(approp_perc)))
# a2<-ds %>% group_by(year) %>%summarise(sd(!is.na(approp_perc)))
# points(ds$year, ds$approp_perc, ylim=c(0.6,1), col="green")
# 
# #pre mean trend through time
# p<-ds %>% group_by(year) %>%summarise(mean(!is.na(pre_perc)))
# p2<-ds %>% group_by(year) %>%summarise(sd(!is.na(pre_perc)))
# plot(p, ylim=c(0.6,1))

#wr_dens trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(wr_dens))
s2<-ds %>% group_by(HUC) %>%summarise(sd(wr_dens))
shp@data<-full_join(shp@data, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","Water_Rights_Density")
spplot(shp[ ,"Water_Rights_Density"])

#rip trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(rip_perc))
s2<-ds %>% group_by(HUC) %>%summarise(sd(rip_perc))
shp@data<-full_join(shp@data, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","Percent_Riparian")
#spplot(shp[ ,"Percent_Riparian"])

#approp trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(approp_perc))
s2<-ds %>% group_by(HUC) %>%summarise(sd(approp_perc))
shp@data<-full_join(shp@data, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC", "Percent_Riparian", "Percent_Post_1914_Appropriative")
#spplot(shp[ ,"Percent_Post-1914_Appropriative"], main="Percent Post-1914 Appropriative")

#pre trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(pre_perc))
s2<-ds %>% group_by(HUC) %>%summarise(sd(pre_perc))
shp@data<-full_join(shp@data, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC", "Percent_Riparian", "Percent_Appropriative", "Percent_Pre_1914")

spplot(shp, c("Percent_Riparian","Percent_Pre_1914", "Percent_Appropriative"), main="Mean Structure of Surface Water Rights")

#FIGURE X
##plot land use for a year at field level
sub<-d
sub<-sub[sub@data$flflag==1,]
sub@data<-sub@data[sub@data$year==14,] #pull out one year at a time
z1<- spplot(sub[, "aglulc"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),
            legendEntries =c(" ","Barren & Fallow", "Grasses", "Grains", "Row Crops","Fruits and Nuts","Uncultivated Cover"),
            main="Land Use Classification in 2014")
z1$legend$inside$args$key$points$cex<-c(0,1,1,1,1,1,1)
plot(z1)

######################
# PLOTS of WR effects on TVP given interactions with time
#gives fixed effects and random slope of fixed effect by HUC 
#this is the same as the intercept of the wr-time interaction effects
# ##################
# stlm3c<-y~ 1  + year + gw + spi + factor(aglulc) + fl_cd+ wr_dens + rip_perc + approp_perc + pre_perc +
#   year*rip_perc +year*approp_perc + year*pre_perc +
#   f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="iid") + f(year1, model="ar1")
# output.stlm3c <- inla (stlm3c,  family = "gaussian", data=ds, verbose=TRUE)#, 
# #control.predictor =list( compute=TRUE))
# saveRDS(output.stlm3c, file='/home/kate/CA/data/output.stlm3c')
# stlm3c<-readRDS('/home/kate/CA/data/output.stlm3c')

stlm3c<-readRDS('output.stlm3c')

b0<-inla.rmarginal(1000,marg=stlm3c$marginals.fixed$`(Intercept)`) #extract the mean effect (intercept) from the fixed marginals
v0k<-matrix(NA,1000,849) #a matrix of the random effect of HUC on intercept (beta0jk)
  for(i in 1:849){
    v0k[ ,i]<-inla.rmarginal(1000,marg=stlm3c$marginals.random$HUC[1:849][[i]])
  }

beta0jk<-b0+v0k #compute beta0jk
beta0jk_quartiles<-t(apply(beta0jk, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
beta0jk_q<-as.data.frame(beta0jk_quartiles)
beta0jk_q<-mutate(beta0jk_q, ID=row.names(beta0jk_q))
beta0jk_q$ID<-as.numeric(beta0jk_q$ID)

#plot variance of the intercept by HUC (Beta0jk)
shp_dup<-shp
shp_dup@data<-full_join(shp_dup@data, beta0jk_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
#spplot(shp_dup, "median", main = "Posterior Median of Random Intercept")#plot the  median intercept
spplot(shp_dup, c("Lower","Median","Upper"), main = "Posterior Median and 95% Credibility Interval of Random Intercept")

#############################################################################
#plot the expected value of the conditional regression of riparian  water rights
#on tvp as a function of time (all other predictors & controls at zero)
#yhat<- (b0 + b2*moderator)+(b1 + b3*moderator)*focal #equation for the conditional regression of focal on y as fcn of moderator
##################################################

#RIPARIAN

moderator <-seq(0,7) # year
focal<-seq(0,1,0.2) #water right percentage
b0<-stlm3c$summary.fixed[1,] #overall intercept
b1<-stlm3c$summary.fixed[12,] #main effect of the focal predictor, rip in this case
b2<-stlm3c$summary.fixed[2,]#main effect of the moderator, year
b3<-stlm3c$summary.fixed[15,] #interaction effect, rip:year

#Plot of the conditional regression of TVP on riparian rights as a fcn of time
plot(focal, ((b0[,1] + b2[,1]*moderator[2])+(b1[,1] + b3[,1]*moderator[2])*focal), type="l", ylim=c(0.42,0.5), xlab="Percent Riparian Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Riparian Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1[,1] + b3[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1[,1] + b3[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[8])+(b1[,1] + b3[,1]*moderator[8])*focal), col="red")
legend(0.7,0.50, c("1", "3", "5", "7"),title ="Years of Drought", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","green","blue", "red"))

 
#plot the confidence in the simple slope of riparian rights on tcp as a function of time
#simple slope <- (b1[,1] + b3[,1]*moderator)*focal[1])
mod<-as.matrix(moderator[c(1,2,4,6,8)])
b1s<-inla.rmarginal(1000,marg=stlm3c$marginals.fixed$rip_perc)
b3s<-inla.rmarginal(1000, marg=stlm3c$marginals.fixed$`year:rip_perc`)
simpleslope<-matrix(NA,1000,5)
simpleslope[,1]<-(b1s + b3s*mod[1])
simpleslope[,2]<-(b1s + b3s*mod[2])
simpleslope[,3]<-(b1s + b3s*mod[3])
simpleslope[,4]<-(b1s + b3s*mod[4])
simpleslope[,5]<-(b1s + b3s*mod[5])
simpleslope_quartiles<-t(apply(simpleslope, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))
plot((mod), simpleslope_quartiles[,2], type="l", ylim=c(-0.01, 0.04), ylab="Simple Slope", xlab="Years of Drought", 
     main="Region of Significance for theSimple Slope of TVP \nRegressed on Percent Riparian Water Rights", cex.main=1)
lines((mod), simpleslope_quartiles[,1], col="red", lty=2)
lines((mod), simpleslope_quartiles[,3], col="red", lty=2)

#APPROPRIATIVE
b0<-stlm3c$summary.fixed[1,] #overall intercept
b1<-stlm3c$summary.fixed[13,] #main effect of the focal predictor, rip in this case
b2<-stlm3c$summary.fixed[2,]#main effect of the moderator, year
b3<-stlm3c$summary.fixed[16,] #interaction effect, rip:year

plot(focal, ((b0[,1] + b2[,1]*moderator[2])+(b1[,1] + b3[,1]*moderator[2])*focal), type="l", ylim=c(0.42,0.51), xlab="Percent Appropriative Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Appropriative Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1[,1] + b3[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1[,1] + b3[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[8])+(b1[,1] + b3[,1]*moderator[8])*focal), col="red")
legend(0.7,0.51, c("1", "3", "5", "7"),title ="Years of Drought", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","green","blue", "red"))

mod<-as.matrix(moderator[c(1,2,4,6,8)])
b1s<-inla.rmarginal(1000,marg=stlm3c$marginals.fixed$approp_perc)
b3s<-inla.rmarginal(1000, marg=stlm3c$marginals.fixed$`year:approp_perc`)
simpleslope<-matrix(NA,1000,5)
simpleslope[,1]<-(b1s + b3s*mod[1])
simpleslope[,2]<-(b1s + b3s*mod[2])
simpleslope[,3]<-(b1s + b3s*mod[3])
simpleslope[,4]<-(b1s + b3s*mod[4])
simpleslope[,5]<-(b1s + b3s*mod[5])
simpleslope_quartiles<-t(apply(simpleslope, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))
plot((mod), simpleslope_quartiles[,2], type="l", ylim=c(-0.02, 0.02), ylab="Simple Slope", xlab="Years of Drought", 
     main="Region of Significance for theSimple Slope of TVP \nRegressed on Percent Appropriative Water Rights", cex.main=1)
lines((mod), simpleslope_quartiles[,1], col="red", lty=2)
lines((mod), simpleslope_quartiles[,3], col="red", lty=2)

#Pre-1914
b0<-stlm3c$summary.fixed[1,] #overall intercept
b1<-stlm3c$summary.fixed[14,] #main effect of the focal predictor, rip in this case
b2<-stlm3c$summary.fixed[2,]#main effect of the moderator, year
b3<-stlm3c$summary.fixed[17,] #interaction effect, rip:year

plot(focal, ((b0[,1] + b2[,1]*moderator[2])+(b1[,1] + b3[,1]*moderator[2])*focal), type="l", ylim=c(0.42,0.54), xlab="Percent Pre-1914 Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Pre-1914 Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1[,1] + b3[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1[,1] + b3[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[8])+(b1[,1] + b3[,1]*moderator[8])*focal), col="red")
legend(0.7,0.54, c("1", "3", "5", "7"),title ="Years of Drought", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","green","blue", "red"))

mod<-as.matrix(moderator[c(1,2,4,6,8)])
b1s<-inla.rmarginal(1000,marg=stlm3c$marginals.fixed$pre_perc)
b3s<-inla.rmarginal(1000, marg=stlm3c$marginals.fixed$`year:pre_perc`)
simpleslope<-matrix(NA,1000,5)
simpleslope[,1]<-(b1s + b3s*mod[1])
simpleslope[,2]<-(b1s + b3s*mod[2])
simpleslope[,3]<-(b1s + b3s*mod[3])
simpleslope[,4]<-(b1s + b3s*mod[4])
simpleslope[,5]<-(b1s + b3s*mod[5])
simpleslope_quartiles<-t(apply(simpleslope, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))
plot((mod), simpleslope_quartiles[,2], type="l", ylim=c(0, 0.1), ylab="Simple Slope", xlab="Years of Drought", 
     main="Region of Significance for theSimple Slope of TVP \nRegressed on Percent Pre-1914 Water Rights", cex.main=1)
lines((mod), simpleslope_quartiles[,1], col="red", lty=2)
lines((mod), simpleslope_quartiles[,3], col="red", lty=2)



#ALL THREE TOGETHER
b0p<-stlm3c$summary.fixed[1,] #overall intercept
b1p<-stlm3c$summary.fixed[14,] #main effect of the focal predictor, pre in this case
b2p<-stlm3c$summary.fixed[2,]#main effect of the moderator, year
b3p<-stlm3c$summary.fixed[17,] #interaction effect, pre:year
b0a<-stlm3c$summary.fixed[1,] #overall intercept
b1a<-stlm3c$summary.fixed[13,] #main effect of the focal predictor, approp in this case
b2a<-stlm3c$summary.fixed[2,]#main effect of the moderator, year
b3a<-stlm3c$summary.fixed[16,] #interaction effect, approp:year
b0r<-stlm3c$summary.fixed[1,] #overall intercept
b1r<-stlm3c$summary.fixed[12,] #main effect of the focal predictor, rip in this case
b2r<-stlm3c$summary.fixed[2,]#main effect of the moderator, year
b3r<-stlm3c$summary.fixed[15,] #interaction effect, rip:year

plot(focal, ((b0p[,1] + b2p[,1]*moderator[1])+(b1p[,1] + b3p[,1]*moderator[1])*focal), type="l", col="red", lwd=2,ylim=c(0.42,0.54), xlab="Percent Pre-1914 Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0p[,1] + b2p[,1]*moderator[4])+(b1p[,1] + b3p[,1]*moderator[4])*focal), col="red", lty=2, lwd=2)
lines(focal, ((b0p[,1] + b2p[,1]*moderator[8])+(b1p[,1] + b3p[,1]*moderator[8])*focal), col="red", lty=4, lwd=2)
lines(focal, ((b0a[,1] + b2a[,1]*moderator[1])+(b1a[,1] + b3a[,1]*moderator[1])*focal), col="black", lwd=2)
lines(focal, ((b0a[,1] + b2a[,1]*moderator[4])+(b1a[,1] + b3a[,1]*moderator[4])*focal), col="black", lty=2, lwd=2)
lines(focal, ((b0a[,1] + b2a[,1]*moderator[8])+(b1a[,1] + b3a[,1]*moderator[8])*focal), col="black", lty=4, lwd=2)
lines(focal, ((b0r[,1] + b2r[,1]*moderator[1])+(b1r[,1] + b3r[,1]*moderator[1])*focal), col="blue", lwd=2)
lines(focal, ((b0r[,1] + b2r[,1]*moderator[4])+(b1r[,1] + b3r[,1]*moderator[4])*focal), col="blue", lty=2, lwd=2)
lines(focal, ((b0r[,1] + b2r[,1]*moderator[8])+(b1r[,1] + b3r[,1]*moderator[8])*focal), col="blue", lty=4, lwd=2)
legend(0.1,0.54, c("Appropriative","Riparian","Pre-1914", "Year 1","Year 4", "Year 8"),title ="Type of Water Right", 
       lwd=c(2,2,2,2,2,2), lty=c(1,1,1,1,2,4), col=c("black","blue", "red","black","black", "black"))

#This plot flips the x axis
plot(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[2]), type="l", col="red",ylim=c(0.43,0.49), xlab="Year of Drought", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[3]), col="red", lty=2)
lines(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[4]), col="red", lty=4)
lines(moderator, ((b0a[,1] + b2a[,1]*moderator)+(b1a[,1] + b3a[,1]*moderator)*focal[2]), col="black")
lines(moderator, ((b0a[,1] + b2a[,1]*moderator)+(b1a[,1] + b3a[,1]*moderator)*focal[3]), col="black", lty=2)
lines(moderator, ((b0a[,1] + b2a[,1]*moderator)+(b1a[,1] + b3a[,1]*moderator)*focal[4]), col="black", lty=4)
lines(moderator, ((b0r[,1] + b2r[,1]*moderator)+(b1r[,1] + b3r[,1]*moderator)*focal[2]), col="blue")
lines(moderator, ((b0r[,1] + b2r[,1]*moderator)+(b1r[,1] + b3r[,1]*moderator)*focal[3]), col="blue", lty=2)
lines(moderator, ((b0r[,1] + b2r[,1]*moderator)+(b1r[,1] + b3r[,1]*moderator)*focal[4]), col="blue", lty=4)
legend(0.1,0.46, c("Appropriative","Riparian","Pre-1914", "20%","40%", "60%"),title ="Type of Water Right", 
       lwd=c(2,2,2,2,2,2), lty=c(1,1,1,1,2,4), col=c("black","blue", "red","black","black", "black"))


###########################################
#Interaction trend plot & map of crossings
############################################

x <- seq(1,8) # Years
#this plot shows us the simple slope of water rights over time - change in effect on tvp over time
#note that this is the average for all spatial units and that the slope for these interactions with time is the same for all HUCs 
plot(x,m1$summary.fixed[16,4]*x + m1$summary.fixed[12,4], type="l", lwd=2, main="",xlab="year",ylab=expression("average effect of water right type on tvp for the central valley"),ylim=c(-0.05,0.05))
lines(x,m1$summary.fixed[17,4]*x + m1$summary.fixed[13,4],lty=1, col="red", lwd=2) #pre1914
lines(x, m1$summary.fixed[15,4]*x + m1$summary.fixed[11,4],lty=1, col="blue", lwd=2) #rip
lines(x,m1$summary.fixed[16,3]*x + m1$summary.fixed[12,3], lty=4, lwd=1) #approp 0.25 confidence 
lines(x,m1$summary.fixed[16,5]*x + m1$summary.fixed[12,5], lty=4, lwd=1) #approp 0.975 confidence
lines(x,m1$summary.fixed[17,3]*x + m1$summary.fixed[13,3],lty=4, col="red", lwd=1) #pre1914 0.25 confidence
lines(x,m1$summary.fixed[17,5]*x + m1$summary.fixed[13,5],lty=4, col="red", lwd=1) #pre1914 0.75 confidence
lines(x, m1$summary.fixed[15,3]*x + m1$summary.fixed[11,3],lty=4, col="blue", lwd=1) #rip 0.25 confidence
lines(x, m1$summary.fixed[15,5]*x + m1$summary.fixed[11,5],lty=4, col="blue", lwd=1) #rip 0.75 confidence
legend(4.5,0.1, c("approp","pre1914"," riparian"), 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","red","blue"))

#join the wr effect info to the huc shapefile
fullplot<-ripplot
fullplot@data<-left_join(fullplot@data, appropplot@data, by="HUC")
fullplot@data<-left_join(fullplot@data, preplot@data, by="HUC")

#calcualte intersection pt (year at which lines for effect of pre and effect of approp cross)
#diff of the intercepts divided by the sum of the slopes
xpt<-fullplot
xpt@data<-mutate(xpt@data, crossings=(xpt@data[,22]-xpt@data[ ,33]) /(m1$summary.fixed[16,4] +m1$summary.fixed[17,4]))
spplot(xpt[ ,35])

######################
# PLOTS of WR effects on TVP given interactions with spi
#gives fixed effects and random slope of fixed effect by HUC 
#this is the same as the intercept of the wr-time interaction effects
##################


#plot variance of the rip_perc slope by HUC, first half (849)of summary.random$HUC is structured, second half is unstructured
rip<- sm1b$summary.random$HUC2[ ,c(1,5)] #pull the median value for rip_perc random effects
rip_struc<-rip[1:849, ] #pull the spatial random effects for rip_perc 
rip_tot<- mutate(rip_struc, fixed_random = (rip_struc[ ,2]  +sm1b$summary.fixed[11,4]), fixed = sm1b$summary.fixed[11,4]) #add both random effects to the median fixed effect
ripplot<-shp
ripplot@data<-full_join(ripplot@data, rip_tot, by=c("HUC"="ID"))
colnames(ripplot@data)<-c("HUC12", "HUC","tvp", "tvp0", "tvp7", "delta_tvp", "rip", "approp", "pre", "random","fixed_random",  "fixed")
spplot(ripplot[ ,10])#plot the  random effects
spplot(ripplot[ ,12])#plot the fixed effects
spplot(ripplot[ ,11])#plot the total effect of rip_perc

#plot variance of the approp_perc slope by HUC, first half (849)of summary.random$HUC is structured, second half is unstructured
approp<- sm1b$summary.random$HUC3[ ,c(1,5)] #pull the median value for approp_perc random effects
approp_struc<-approp[1:849, ] #pull the first set of random effects for approp_perc (structured effects)
approp_tot<- mutate(approp_struc, fixed_random = (approp_struc[ ,2]  +sm1b$summary.fixed[12,4]),  fixed= sm1b$summary.fixed[12,4]) #add both random effects to the median fixed effect
appropplot<-shp
appropplot@data<-full_join(appropplot@data, approp_tot, by=c("HUC"="ID"))
colnames(appropplot@data)<-c("HUC12", "HUC","tvp", "tvp0", "tvp7", "delta_tvp", "rip", "approp", "pre", "randoma","fixed_randoma",  "fixeda")
spplot(appropplot[ ,10])#plot the  random effects
spplot(appropplot[ ,12])#plot the fixedeffects
spplot(appropplot[ ,11])#plot the total effect of approp_perc

#plot variance of the pre_perc slope by HUC, first half (849)of summary.random$HUC is structured, second half is unstructured
pre<- sm1b$summary.random$HUC4[ ,c(1,5)] #pull the median value for pre_perc random effects
pre_struc<-pre[1:849, ] #pull the first set of random effects for pre_perc (structured effects)
pre_tot<- mutate(pre_struc, fixed_random = (pre_struc[ ,2]  +sm1b$summary.fixed[13,4]), fixed = sm1b$summary.fixed[13,4]) #add both random effects to the median fixed effect
preplot<-shp
preplot@data<-full_join(preplot@data, pre_tot, by=c("HUC"="ID"))
colnames(preplot@data)<-c("HUC12", "HUC","tvp", "tvp0", "tvp7", "delta_tvp", "rip", "approp", "pre", "randomp","fixed_randomp","fixedp")
spplot(preplot[ ,10])#plot the structured random effects
spplot(preplot[ ,12])#plot the unstructured random effects
spplot(preplot[ ,11])#plot the total effect of pre_perc

####################
#same interaction plot for spi
x <- seq(-5,5) # spi values
plot(x,sm1b$summary.fixed[16,1]*x + sm1b$summary.fixed[12,1], type="l", lwd=2, main="",xlab="SPI",ylab=expression("average effect of water right type on tvp for the central valley"),ylim=c(-0.25,0.25))
lines(x,sm1b$summary.fixed[17,1]*x + sm1b$summary.fixed[13,1],lty=1, col="red", lwd=2) #pre1914
lines(x, sm1b$summary.fixed[15,1]*x + sm1b$summary.fixed[11,1],lty=1, col="blue", lwd=2) #rip
lines(x,sm1b$summary.fixed[16,3]*x + sm1b$summary.fixed[12,3], lty=4, lwd=1) #approp 0.25 confidence 
lines(x,sm1b$summary.fixed[16,5]*x + sm1b$summary.fixed[12,5], lty=4, lwd=1) #approp 0.975 confidence
lines(x,sm1b$summary.fixed[17,3]*x + sm1b$summary.fixed[13,3],lty=4, col="red", lwd=1) #pre1914 0.25 confidence
lines(x,sm1b$summary.fixed[17,5]*x + sm1b$summary.fixed[13,5],lty=4, col="red", lwd=1) #pre1914 0.75 confidence
lines(x, sm1b$summary.fixed[15,3]*x + sm1b$summary.fixed[11,3],lty=4, col="blue", lwd=1) #rip 0.25 confidence
lines(x, sm1b$summary.fixed[15,5]*x + sm1b$summary.fixed[11,5],lty=4, col="blue", lwd=1) #rip 0.75 confidence
legend(2,0.5, c("approp","pre1914"," riparian"), 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","red","blue"))

#join the wr effect info to the huc shapefile
fullplot<-ripplot
fullplot@data<-left_join(fullplot@data, appropplot@data, by="HUC")
fullplot@data<-left_join(fullplot@data, preplot@data, by="HUC")

#calcualte intersection pt (year at which lines for effect of pre and effect of approp cross)
#diff of the intercepts divided by the sum of the slopes
xpt<-fullplot
xpt@data<-mutate(xpt@data, crossings=(xpt@data[,22]-xpt@data[ ,33]) /(sm1b$summary.fixed[16,4] +sm1b$summary.fixed[17,4]))
spplot(xpt[ ,35])
#crossing for approp and rip
xpt@data<-mutate(xpt@data, crossings2=(xpt@data[,22]-xpt@data[ ,11]) /(sm1b$summary.fixed[16,4] +sm1b$summary.fixed[15,4]))
spplot(xpt[ ,36])

######################################
#Plot the time trend
######################################
x <- seq(1,8) # Years
plot(x,m1$summary.fixed[14,4]*x, type="l", main="",xlab="t",ylab=expression(beta*t), ylim=c(-0.15,0.15))
lines(m1$summary.fixed[14,3]*x,lty=2)
lines(m1$summary.fixed[14,5]*x,lty=2)
lines(m1$summary.random$year1[ ,1]+1, m1$summary.random$year1[ ,5], col="red")
lines(m1$summary.random$year1[ ,1]+1, m1$summary.random$year1[ ,6], col="red", lty=2)
lines(m1$summary.random$year1[ ,1]+1, m1$summary.random$year1[ ,4], col="red", lty=2)

#############################################################
#Plot the space trend
###############################################
spacetrend<-shp
random<- m1d2$summary.random$HUC[ ,c(1,5)] #pull the median bym random effects
struc<-random[1:849, ] #pull the spatially structured random effects  
unstruc<-random[850:1698,] #pull the unstructured spatial effects
totrandom<- mutate(struc, tot_random = (struc[ ,2]  +unstruc[ ,2]), unstruc=unstruc[ ,2]) #add both random effects to the median fixed effect
spacetrend@data<-full_join(spacetrend@data, totrandom, by=c("HUC"="ID"))
colnames(spacetrend@data)<-c("HUC12", "HUC","tvp", "tvp0", "tvp7", "delta_tvp", "rip", "approp", "pre","structured" ,"tot_random", "unstructured" )
spplot(spacetrend[ ,10])#plot the  spatially structured random effects
spplot(spacetrend[ ,12])#plot the  spatially  unstructured random effects
spplot(spacetrend[ ,11])#plot the  total spatial random effects

#####################################
#Plot the standardized residuals

res<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
res<-cbind(res,stlm3cfv$summary.fitted.values$mean,stlm3cfv$summary.fitted.values$sd)
colnames(res)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
res<-mutate(res, residuals = (px_tvp-mean)/sd)
coordinates(res)<- ~lon+lat
res@proj4string <- CRS('+init=epsg:3310')

res_dup<-res
res@data<-res@data[res@data$year==0,]
p1<-spplot(res[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="Standardized \nResiduals: 2007", cuts=c(-155,-34,-17,17,34,155))
p1$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
plot(p1)

res<-res_dup
res@data<-res@data[res@data$year==7,]
p2<-spplot(res[ ,"residuals"],cex=0.1,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="Standardized \nResiduals: 2014",cuts=c(-155,-34,-17,17,34,155))
p2$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
plot(p2)

res<-res_dup
res@data<-res@data[res@data$year==4,]
p3<-spplot(res[ ,"residuals"],cex=0.1,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="Standardized \n Residuals: 2010",cuts=c(-155,-34,-17,17,34,155))
p3$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
plot(p3)

print(p1, position=c(0,0,0.5,1), more=TRUE)
print(p2, position=c(0.5,0,1,1))

##################
#residuals through time
##############
r<-res@data
t<-r %>% group_by(year) %>% summarise(mean=mean(residuals))
t2<-r %>% group_by(year) %>%summarise(sd=sd(residuals))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim = c(-18,18), main="Mean Residuals through Time", xlab="Years since 2007", ylab="Mean Residuals for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)


##################
#fitted values (estimates of tvp) through time
##############
r<-res@data
t<-r %>% group_by(year) %>% summarise(mean=mean(mean))
t2<-r %>% group_by(year) %>%summarise(sd=mean(sd))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim = c(0,1), main="Mean Fitted Values through Time", xlab="Years since 2007", ylab="Mean Fitted Values for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

#####################################
#Plot the fitted values
####################################

res<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
res<-cbind(res,stlm3cfv$summary.fitted.values$mean,stlm3cfv$summary.fitted.values$sd)
colnames(res)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
res<-mutate(res, residuals = (px_tvp-mean)/sd)
coordinates(res)<- ~lon+lat
res@proj4string <- CRS('+init=epsg:3310')

res_dup<-res
res@data<-res@data[res@data$year==0,]
p1<-spplot(res[ ,"mean"],cex=0.2,key.space= list(x=0.4, y=0.95, corner=c(0,1)), main="Fitted Values: 2007",cuts=c(0,0.31,0.45,0.6,0.73,1))
p1$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
plot(p1)

res<-res_dup
res@data<-res@data[res@data$year==7,]
p2<-spplot(res[ ,"mean"],cex=0.1,key.space= list(x=0.4, y=0.95, corner=c(0,1)), main="Fitted Values: 2014",cuts=c(0,0.31,0.45,0.6,0.73,1))
p2$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
plot(p2)

res<-res_dup
res@data<-res@data[res@data$year==4,]
p3<-spplot(res[ ,"mean"],cex=0.1,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="Fitted Values: 2010")#,cuts=c(0,0.31,0.43,0.67,0.79,1))
p3$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
plot(p3)

print(p1, position=c(0,0,0.5,1), more=TRUE)
print(p2, position=c(0.5,0,1,1))
