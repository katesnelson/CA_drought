library(memisc)
library(dplyr)
library (gstat)
library(sp)
library(nlme)
library(rgdal)
library(INLA)
library(spdep)
library(raster)


# setwd('C:/Users/nelsonks/Desktop/drought')
# setwd('E:/drought')
setwd('/home/kate/R/x86_64-pc-linux-gnu-library/3.3/INLA/')
##############################################
#Read in datasets to use for inla models
##############################################
#full cleaned pixel dataset of central valley as spatial points data frame
#d <- readRDS('/data/emily/WF/kate/final_data/cv.dat.clean.rds')

###  full dataset with lat lon columns and updated wr info --> check for cleanliness 
d<- readRDS('/home/kate/CA/data/cv.dat2.rds') 
#d<- readRDS('cv.dat2.rds') 
#plot(d)


#full central valley watershed as spatial polygons data frame
#cv <- shapefile('/data/emily/WF/kate/shp/cv_huc.shp')
cv <- shapefile('/home/kate/CA/data/cv_huc.shp')
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
#temp <- poly2nb(shp, queen=T)

#convert to a sparse matrix to reduce memory
H.adj <- nb2mat(temp, style ="B") #convert neighbor list to binary coded neighbor weights matrix
H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse sytle matrix

#plot spatial object
image(inla.graph2matrix(H.adj), xlab="", ylab="")

# #write out list as an INLA object to wd
# nb2INLA("/home/kate/CA/data/H.adj", H.adj)

# #############################################
# #Create SPDE object
# ##############################################
# coords <- ds[ ,c(4,3)] #create a list of coordinates from lat,lon in ds
# coords<-unique(coords) # use unique coords to create mesh
# mesh<-inla.mesh.2d(loc=coords,  max.edge=200000, cutoff = 5000 ) #create a mesh for estimation, trouble runing when specify an offset
# plot(mesh) #check projection stuff
# spde <-inla.spde2.matern(mesh=mesh, alpha=2) #pg 206, this mesh is for accounting for continuous spatial processes (landscape,etc..), and not for interpolating predictor values
# # creating the observation/projection matrix (sparse weight matrix pg 205)
# 
# coords.alltime<-as.matrix(coords[ds$POINTID, c("lon", "lat")]) #list of coordinates for all times
# #y<-ds$year + 1
# #n.year<-length(unique(y))
# A.est <-inla.spde.make.A(mesh=mesh,loc=coords.alltime)#, group=y, n.group=n.year) #object that reconciles mesh and observations
# 
# s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)

########################################
#Prep data for analysis
#######################################

y=ds$px_tvp #dependent variable

ds$year <- ds$year - 7 #center time at 2007 = 0

#build duplicate indexes for model building
ds <- mutate(ds, year1=year, year2=year, year3=year, year4 = year, year5 = year, year6 = year) #add duplicate year grouping index
ds<-mutate(ds,POINTID1=POINTID)#add duplicate pixel grouping index
ds <- mutate(ds, HUC1=HUC,HUC2=HUC, HUC3=HUC, HUC4=HUC, HUC5=HUC, HUC6=HUC, HUC7=HUC, HUC8=HUC,HUC9=HUC, 
             HUC10=HUC, HUC11=HUC, HUC12=HUC, HUC13=HUC, HUC14=HUC, HUC15=HUC, HUC16=HUC, HUC17=HUC) #add duplicates of HUC index for random slope grouping
ds<- mutate(ds, wr_dens1=wr_dens, rip_perc1=rip_perc, approp_perc1=approp_perc, pre_perc1=pre_perc)
ds$agluid<-as.numeric(ds$aglulc)#landuse id var for function fields
#s.index1 <- inla.spde.make.index(name="spatial.field1", n.spde=spde$n.spde)

#save a copy pre-standardization for later use
dsdup<-ds

#standardize independent vars to ease interpretation

ds$gw<- (ds$gw - mean(ds$gw, na.rm=T))/sd(ds$gw, na.rm=T)
ds$spi<- (ds$spi - mean(ds$spi, na.rm=T))/sd(ds$spi, na.rm=T) 
ds$fl_cd<-(ds$fl_cd - mean(ds$fl_cd, na.rm=TRUE))/sd(ds$fl_cd, na.rm=TRUE)
ds$wr_dens<- (ds$wr_dens - mean(ds$wr_dens, na.rm=T))/sd(ds$wr_dens,na.rm=T)                             
ds$pre_perc<-ds$pre_perc*100
ds$rip_perc<-ds$rip_perc*100
ds$approp_perc<-ds$approp_perc*100
ds$pre_perc<- (ds$pre_perc - mean(ds$pre_perc, na.rm=T))/sd(ds$pre_perc, na.rm=T)
ds$rip_perc<- (ds$rip_perc - mean(ds$rip_perc, na.rm=T))/sd(ds$rip_perc, na.rm=T)
ds$approp_perc<- (ds$approp_perc - mean(ds$approp_perc, na.rm=T))/sd(ds$approp_perc, na.rm=T)
ds$px_tvp<- (ds$px_tvp - mean(ds$px_tvp, na.rm=T))/sd(ds$px_tvp, na.rm=T)

y<-ds$px_tvp

#dataset broken out by landuse category
ds_E1<-dsdup[dsdup$aglulc==1, ]

ds_E1$gw<- (ds_E1$gw - mean(ds_E1$gw, na.rm=T))/sd(ds_E1$gw, na.rm=T)
ds_E1$spi<- (ds_E1$spi - mean(ds_E1$spi, na.rm=T))/sd(ds_E1$spi, na.rm=T) 
ds_E1$fl_cd<-(ds_E1$fl_cd - mean(ds_E1$fl_cd, na.rm=TRUE))/sd(ds_E1$fl_cd, na.rm=TRUE)
ds_E1$wr_dens<- (ds_E1$wr_dens - mean(ds_E1$wr_dens, na.rm=T))/sd(ds_E1$wr_dens,na.rm=T)                             
ds_E1$pre_perc<-ds_E1$pre_perc*100
ds_E1$rip_perc<-ds_E1$rip_perc*100
ds_E1$approp_perc<-ds_E1$approp_perc*100
ds_E1$pre_perc<- (ds_E1$pre_perc - mean(ds_E1$pre_perc, na.rm=T))/sd(ds_E1$pre_perc, na.rm=T)
ds_E1$rip_perc<- (ds_E1$rip_perc - mean(ds_E1$rip_perc, na.rm=T))/sd(ds_E1$rip_perc, na.rm=T)
ds_E1$approp_perc<- (ds_E1$approp_perc - mean(ds_E1$approp_perc, na.rm=T))/sd(ds_E1$approp_perc, na.rm=T)
ds_E1$px_tvp<- (ds_E1$px_tvp - mean(ds_E1$px_tvp, na.rm=T))/sd(ds_E1$px_tvp, na.rm=T)
y_E1<-ds_E1$px_tvp

#
ds_E2<-dsdup[dsdup$aglulc==2, ]

ds_E2$gw<- (ds_E2$gw - mean(ds_E2$gw, na.rm=T))/sd(ds_E2$gw, na.rm=T)
ds_E2$spi<- (ds_E2$spi - mean(ds_E2$spi, na.rm=T))/sd(ds_E2$spi, na.rm=T) 
ds_E2$fl_cd<-(ds_E2$fl_cd - mean(ds_E2$fl_cd, na.rm=TRUE))/sd(ds_E2$fl_cd, na.rm=TRUE)
ds_E2$wr_dens<- (ds_E2$wr_dens - mean(ds_E2$wr_dens, na.rm=T))/sd(ds_E2$wr_dens,na.rm=T)                             
ds_E2$pre_perc<-ds_E2$pre_perc*100
ds_E2$rip_perc<-ds_E2$rip_perc*100
ds_E2$approp_perc<-ds_E2$approp_perc*100
ds_E2$pre_perc<- (ds_E2$pre_perc - mean(ds_E2$pre_perc, na.rm=T))/sd(ds_E2$pre_perc, na.rm=T)
ds_E2$rip_perc<- (ds_E2$rip_perc - mean(ds_E2$rip_perc, na.rm=T))/sd(ds_E2$rip_perc, na.rm=T)
ds_E2$approp_perc<- (ds_E2$approp_perc - mean(ds_E2$approp_perc, na.rm=T))/sd(ds_E2$approp_perc, na.rm=T)
ds_E2$px_tvp<- (ds_E2$px_tvp - mean(ds_E2$px_tvp, na.rm=T))/sd(ds_E2$px_tvp, na.rm=T)
y_E2<-ds_E2$px_tvp

#
ds_E3<-dsdup[dsdup$aglulc==3, ]

ds_E3$gw<- (ds_E3$gw - mean(ds_E3$gw, na.rm=T))/sd(ds_E3$gw, na.rm=T)
ds_E3$spi<- (ds_E3$spi - mean(ds_E3$spi, na.rm=T))/sd(ds_E3$spi, na.rm=T) 
ds_E3$fl_cd<-(ds_E3$fl_cd - mean(ds_E3$fl_cd, na.rm=TRUE))/sd(ds_E3$fl_cd, na.rm=TRUE)
ds_E3$wr_dens<- (ds_E3$wr_dens - mean(ds_E3$wr_dens, na.rm=T))/sd(ds_E3$wr_dens,na.rm=T)                             
ds_E3$pre_perc<-ds_E3$pre_perc*100
ds_E3$rip_perc<-ds_E3$rip_perc*100
ds_E3$approp_perc<-ds_E3$approp_perc*100
ds_E3$pre_perc<- (ds_E3$pre_perc - mean(ds_E3$pre_perc, na.rm=T))/sd(ds_E3$pre_perc, na.rm=T)
ds_E3$rip_perc<- (ds_E3$rip_perc - mean(ds_E3$rip_perc, na.rm=T))/sd(ds_E3$rip_perc, na.rm=T)
ds_E3$approp_perc<- (ds_E3$approp_perc - mean(ds_E3$approp_perc, na.rm=T))/sd(ds_E3$approp_perc, na.rm=T)
ds_E3$px_tvp<- (ds_E3$px_tvp - mean(ds_E3$px_tvp, na.rm=T))/sd(ds_E3$px_tvp, na.rm=T)
y_E3<-ds_E3$px_tvp

#
ds_E4<-dsdup[dsdup$aglulc==4, ]

ds_E4$gw<- (ds_E4$gw - mean(ds_E4$gw, na.rm=T))/sd(ds_E4$gw, na.rm=T)
ds_E4$spi<- (ds_E4$spi - mean(ds_E4$spi, na.rm=T))/sd(ds_E4$spi, na.rm=T) 
ds_E4$fl_cd<-(ds_E4$fl_cd - mean(ds_E4$fl_cd, na.rm=TRUE))/sd(ds_E4$fl_cd, na.rm=TRUE)
ds_E4$wr_dens<- (ds_E4$wr_dens - mean(ds_E4$wr_dens, na.rm=T))/sd(ds_E4$wr_dens,na.rm=T)                             
ds_E4$pre_perc<-ds_E4$pre_perc*100
ds_E4$rip_perc<-ds_E4$rip_perc*100
ds_E4$approp_perc<-ds_E4$approp_perc*100
ds_E4$pre_perc<- (ds_E4$pre_perc - mean(ds_E4$pre_perc, na.rm=T))/sd(ds_E4$pre_perc, na.rm=T)
ds_E4$rip_perc<- (ds_E4$rip_perc - mean(ds_E4$rip_perc, na.rm=T))/sd(ds_E4$rip_perc, na.rm=T)
ds_E4$approp_perc<- (ds_E4$approp_perc - mean(ds_E4$approp_perc, na.rm=T))/sd(ds_E4$approp_perc, na.rm=T)
ds_E4$px_tvp<- (ds_E4$px_tvp - mean(ds_E4$px_tvp, na.rm=T))/sd(ds_E4$px_tvp, na.rm=T)
y_E4<-ds_E4$px_tvp

#
ds_E5<-dsdup[dsdup$aglulc==5, ]

ds_E5$gw<- (ds_E5$gw - mean(ds_E5$gw, na.rm=T))/sd(ds_E5$gw, na.rm=T)
ds_E5$spi<- (ds_E5$spi - mean(ds_E5$spi, na.rm=T))/sd(ds_E5$spi, na.rm=T) 
ds_E5$fl_cd<-(ds_E5$fl_cd - mean(ds_E5$fl_cd, na.rm=TRUE))/sd(ds_E5$fl_cd, na.rm=TRUE)
ds_E5$wr_dens<- (ds_E5$wr_dens - mean(ds_E5$wr_dens, na.rm=T))/sd(ds_E5$wr_dens,na.rm=T)                             
ds_E5$pre_perc<-ds_E5$pre_perc*100
ds_E5$rip_perc<-ds_E5$rip_perc*100
ds_E5$approp_perc<-ds_E5$approp_perc*100
ds_E5$pre_perc<- (ds_E5$pre_perc - mean(ds_E5$pre_perc, na.rm=T))/sd(ds_E5$pre_perc, na.rm=T)
ds_E5$rip_perc<- (ds_E5$rip_perc - mean(ds_E5$rip_perc, na.rm=T))/sd(ds_E5$rip_perc, na.rm=T)
ds_E5$approp_perc<- (ds_E5$approp_perc - mean(ds_E5$approp_perc, na.rm=T))/sd(ds_E5$approp_perc, na.rm=T)
ds_E5$px_tvp<- (ds_E5$px_tvp - mean(ds_E5$px_tvp, na.rm=T))/sd(ds_E5$px_tvp, na.rm=T)
y_E5<-ds_E5$px_tvp

#
ds_E6<-dsdup[dsdup$aglulc==6, ]

ds_E6$gw<- (ds_E6$gw - mean(ds_E6$gw, na.rm=T))/sd(ds_E6$gw, na.rm=T)
ds_E6$spi<- (ds_E6$spi - mean(ds_E6$spi, na.rm=T))/sd(ds_E6$spi, na.rm=T) 
ds_E6$fl_cd<-(ds_E6$fl_cd - mean(ds_E6$fl_cd, na.rm=TRUE))/sd(ds_E6$fl_cd, na.rm=TRUE)
ds_E6$wr_dens<- (ds_E6$wr_dens - mean(ds_E6$wr_dens, na.rm=T))/sd(ds_E6$wr_dens,na.rm=T)                             
ds_E6$pre_perc<-ds_E6$pre_perc*100
ds_E6$rip_perc<-ds_E6$rip_perc*100
ds_E6$approp_perc<-ds_E6$approp_perc*100
ds_E6$pre_perc<- (ds_E6$pre_perc - mean(ds_E6$pre_perc, na.rm=T))/sd(ds_E6$pre_perc, na.rm=T)
ds_E6$rip_perc<- (ds_E6$rip_perc - mean(ds_E6$rip_perc, na.rm=T))/sd(ds_E6$rip_perc, na.rm=T)
ds_E6$approp_perc<- (ds_E6$approp_perc - mean(ds_E6$approp_perc, na.rm=T))/sd(ds_E6$approp_perc, na.rm=T)
ds_E6$px_tvp<- (ds_E6$px_tvp - mean(ds_E6$px_tvp, na.rm=T))/sd(ds_E6$px_tvp, na.rm=T)
y_E6<-ds_E6$px_tvp

#Dependent var for Bernoulli likelihood, logistic model 
y_H<-ds$aglulc
y_H[y_H!=1 ]<-0
y_H<-as.character(y_H)
y_H<-as.data.frame(as.numeric(y_H)) #barren and fallow is 1, all other landuse cats are 0
##################################################
#Models 
##################################################


##Null 3 level MLM ## 
  #yijk = b000 + u00k + u0jk + eijk 
NULL_MOD <-y ~ 1 + f(HUC, model="iid") + f(POINTID, model = "iid")
NULL_OUT <- inla (NULL_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(NULL_OUT, file='/home/kate/CA/data/N')
N<-readRDS('/home/kate/CA/data/N')
#DIC is  699735.95

##Random int 3-level MLM ## 
  #yijk = b000 + u00k + u0jk + eijk, where u00k = v00k + s00k (area specific and spatially structured randomness)
A_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid")
A_OUT <- inla (A_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(A_OUT, file='/home/kate/CA/data/A')
A<-readRDS('/home/kate/CA/data/A')
  #DIC (699749.79) is not much different from model N
  #15 hr runtime

##Random int 3-level MLM with correlated residuals ## 
#yijk = b000 + u00k + u0jk + eijk, where u00k = v00k + s00k (area specific and spatially structured randomness)
A3_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + f(POINTID1, HUC1, model = "iid")
A3_OUT <- inla (A3_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(A3_OUT, file='/home/kate/CA/data/A3')
A3<-readRDS('/home/kate/CA/data/A3')
#DIC (699890.39) is not improved from model A
#15 hr runtime

##Random int 3-level ecological model ## 
  #yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + eijk, where u00k = v00k + s00k
B_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") +
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc)
B_OUT <- inla (B_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/B')
B<-readRDS('/home/kate/CA/data/B')
  #DIC (678037.67) is improved (decreased) over model A
  #16.5 hr runtime

##Random int 3-level ecological model with predictor spatial effects (random slope) ## 
#yijk = b000 + v00k + s00k + u0jk + b10k + s10k + b20k + s20k + b30k + s30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + eijk
C_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + #f(POINTID, model = "iid") +
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) + 
  f(HUC1,rip_perc, model="besag", graph=H.adj, scale.model=T) +
  f(HUC2,pre_perc, model="besag", graph=H.adj, scale.model=T)+ 
  f(HUC3,approp_perc, model="besag", graph=H.adj, scale.model=T)
C_OUT <- inla (C_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE)) #crashes 
saveRDS(C_OUT, file='/home/kate/CA/data/C')
C<-readRDS('/home/kate/CA/data/C')
#crashes, regardless of modifications

##Random int 3-level ecological model with spi-wr interactions## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
D_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") +
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
D_OUT <- inla (D_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(D_OUT, file='/home/kate/CA/data/D')
D<-readRDS('/home/kate/CA/data/D')
#DIC (677916.30) is improved (decreased) slightly over model B
#17+ hr runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
E_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
E_OUT <- inla (E_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(E_OUT, file='/home/kate/CA/data/E')
E<-readRDS('/home/kate/CA/data/E')
#DIC (538590.29) is much improved (decreased) over model D
#18+ hr runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect and random temporal effect ## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
F_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + f(year, model = "rw1") +
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
F_OUT <- inla (F_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(F_OUT, file='/home/kate/CA/data/F')
ModelF<-readRDS('/home/kate/CA/data/F')
#DIC (542463.61) is NOT improved (increase) over model E
#30+ hr runtime
#no reason to add a random time effect in addition to a fixed effect

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
E_MOD2 <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") +  f(year, model = "rw1") +
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
E_OUT2 <- inla (E_MOD2, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(E_OUT2, file='/home/kate/CA/data/E2')
E2<-readRDS('/home/kate/CA/data/E2')
#DIC (538971.39) is not much different from model E
#coefficients are about the same as model E so no substantive reason for choosing model with random time effect over one with fixed time effect

##Random int 3-level ecological model with spi-wr-landuse interactions and year fixed-effect## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) + factor(year) +
  rip_perc*spi*factor(aglulc) + pre_perc*spi*factor(aglulc) + approp_perc*spi*factor(aglulc)
G_OUT <- inla (G_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT, file='/home/kate/CA/data/G')
G<-readRDS('/home/kate/CA/data/G')
#DIC (537409.1) is  improved (decreased) over model E
#20.8+ hr runtime
#Practically non-interpretable

##Non-spatial, non-nested model with spi-wr interactions and year fixed-effect check## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
A2 <-y ~ 1 + factor(year) +
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(aglulc) + 
  rip_perc*spi + pre_perc*spi + approp_perc*spi
A2 <- inla (A2, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(A2, file='/home/kate/CA/data/A2')
A2<-readRDS('/home/kate/CA/data/A2')
#DIC (1252412.78) is not improved (increased) over model E or model A 
#spacetime trends only better explain data than linear predictors only
# and spacetime plus linear predictors better explains data than just spacetime, just linear predictors, and other combinations of predictors
#note that size of wr predictor effects is larger than in Model E, indicating that some of the effect of water rights may be due simply to location
#or...that the spatial effect soaks up some of the water right effect
#by including the spatial effect we are trying to isolate the part of the effect that is due to the legal structure and not the location of water rights

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect for land-use 1## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD1 <-y_E1 ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
G_OUT1 <- inla (G_MOD1, family = "gaussian", data=ds_E1, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT1, file='/home/kate/CA/data/G1')
G1<-readRDS('/home/kate/CA/data/G1')
#DIC (43416.00)
#6 min runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect for land-use 2## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD2 <-y_E2 ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
G_OUT2 <- inla (G_MOD2, family = "gaussian", data=ds_E2, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT2, file='/home/kate/CA/data/G2')
G2<-readRDS('/home/kate/CA/data/G2')
#DIC (52566.85)
#6 min runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect for land-use 3## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD3 <-y_E3 ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
G_OUT3 <- inla (G_MOD3, family = "gaussian", data=ds_E3, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT3, file='/home/kate/CA/data/G3')
G3<-readRDS('/home/kate/CA/data/G3')
#DIC (76876.80)
#12 min runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect for land-use 4## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD4 <-y_E4 ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
G_OUT4 <- inla (G_MOD4, family = "gaussian", data=ds_E4, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT4, file='/home/kate/CA/data/G4')
G4<-readRDS('/home/kate/CA/data/G4')
#DIC (45410.69)
#12 min runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect for land-use 5## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD5 <-y_E5 ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
G_OUT5 <- inla (G_MOD5, family = "gaussian", data=ds_E5, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT5, file='/home/kate/CA/data/G5')
G5<-readRDS('/home/kate/CA/data/G5')
#DIC (110656.17)
#23 min runtime

##Random int 3-level ecological model with spi-wr interactions and year fixed-effect for land-use 6## 
#yijk = b000 + u00k + u0jk + beta10k + beta20k + beta30k + beta40k + beta5jk + beta60k + beta70k + beta8jk + beta90k + beta100k + beta110k+ eijk, where u00k = v00k + s00k
G_MOD6 <-y_E6 ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
G_OUT6 <- inla (G_MOD6, family = "gaussian", data=ds_E6, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(G_OUT6, file='/home/kate/CA/data/G6')
G6<-readRDS('/home/kate/CA/data/G6')
#DIC (244347.41)
#5.4 hr runtime

#Random int 3-level ecological model with spi-wr interactions and year fixed-effect for Bernoulli likelihood of landuse cat being barren&fallow
H_MOD <-y_H ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE) + f(POINTID, model = "iid") + 
  rip_perc + pre_perc + approp_perc + wr_dens + gw + spi + fl_cd + factor(year) +
  rip_perc*spi + pre_perc*spi + approp_perc*spi
H_OUT <- inla (H_MOD, family = "binomial", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(H_OUT, file='/home/kate/CA/data/H')
H<-readRDS('/home/kate/CA/data/H')
#DIC(156576.57)
#12.6 hrs
#see page 142 for exponentiating for interp. use antilogit for intercept and exponentiate for betas
# inla.emarginal on coeffs indicates that when spi = 0 more wr reduces prob of b& f abd pre reduces it most, 
#increases in spi reduce prob of b&f, gw recharge increases prob of b&f slightly, and 
#as spi increases (less severe drought) more rip and pre further decrease prob of b&f and approp increases prob of b&f ?
#note that yr fixed effect increases prob the most

#check the size of space effect and spi effect when adding time fe
###################################
#Null Model Diagnostics
###################################
N<-readRDS('/home/kate/CA/data/N')

#extract precisions for random effects, convert to variance, and extract expected value (ie mean variance)
prec<-N$marginals.hyperpar$`Precision for the Gaussian observations`
eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-N$marginals.hyperpar$`Precision for HUC` 
u00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2) 
prec3<-N$marginals.hyperpar$`Precision for POINTID` 
u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)

#calculate interclass correlation coefficients
ICC1<-(eijk/(eijk+u0jk+u00k))
ICC2<-(u0jk/(u0jk+u00k+eijk))
ICC3<-(u00k/(u00k+u0jk+eijk))

#compute random intercept
  #b0jk = b000 + u00k + u0jk

b000<-inla.rmarginal(1000,marg=N$marginals.fixed$`(Intercept)`) #random draw for the marginal distribution of the mean effect (intercept) 

#setup id sequence needed for random draw loop
id2<-as.data.frame(unique(ds$POINTID))#list of unique POINTIDs (level 2 units)
n2<- seq(1,length(id2$`unique(ds$POINTID)`))#count of level 2 units
id2<-cbind(n2,id2) 
ds<-left_join(ds,id2, by=c("POINTID"= "unique(ds$POINTID)")) #create a column in the dataset that links the original POINTID to the simplified sequence

#build a matrix of level 2 and level 3 combined random effects
u0jk<-matrix(NA,1000,length(n2))#an empty matrix to hold the random effects of POINTID and HUC on the intercept term b0jk
U<-matrix(NA,1000,length(n2)) #empty matrix to hold the level 2 and level 3 random effects
for(i in 1:length(n2)){
   u0jk[ ,i]<-inla.rmarginal(1000,marg=N$marginals.random$POINTID[[i]]) #random draw from the marginal distribtuion of the level 2 RE
   n3=first(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located                       
   u00k<-inla.rmarginal(1000,marg=N$marginals.random$HUC[[n3]]) #random draw from the marginal distribution of the level 3 RE for the group in which the level 2 unit is nested
   U[ ,i]<- u0jk[ ,i] + u00k #add the appropriate level 3 random effects to the level 2 effects
}
  
beta0jk<-b000 + U #compute beta0jk
beta0jk_quartiles<-t(apply(beta0jk, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
beta0jk_q<-as.data.frame(beta0jk_quartiles)
beta0jk_q<-mutate(beta0jk_q, ID=row.names(beta0jk_q))
beta0jk_q$ID<-as.numeric(beta0jk_q$ID)

#plot the random intercept (beta0jk)
ds_dup<-ds
ds_dup<-left_join(ds, beta0jk_q, by=c("n2"= "ID"))
coordinates(ds_dup)= c("lon", "lat")
Nplot<-ds_dup[ ,74:77]
colnames(Nplot@data)<-c("n2","Lower", "Median", "Upper")
Nplot@proj4string <- CRS('+init=epsg:3310')
spplot(Nplot, c("Lower","Median","Upper"), main = "N: Posterior Median and 95% Credibility Interval of Random Intercept", cex =0.1)
plotNull<- spplot(Nplot, c("Median"), key.space= list(x=0.5, y=0.95, corner=c(0,1)), main = "N: Posterior Median of Random Intercept", cex =0.2)
plotNull$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
print(plotNull)
  #random intercept shows spatial clustering on an ecological scale (i.e. valley vs mountain)

#plot the standardized residuals --> note that a positive residual indicates that the mean fitted value underestimates the observed value
Nres<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
Nres<-cbind(Nres,N$summary.fitted.values$mean,N$summary.fitted.values$sd)
colnames(Nres)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
Nres<-mutate(Nres, residuals = (px_tvp-mean)/sd)
plot(Nres$year,Nres$residuals, cex=0.5) #plot standardized residuals through time
  #standardized residuals have large spread in a single year and a smaller sinusoidal variance across time

coordinates(Nres)<- ~lon+lat
Nres@proj4string <- CRS('+init=epsg:3310')
Nres@data<-Nres@data[Nres@data$year==0,]
plotNullres<-spplot(Nres[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="N:Standardized \nResiduals: 2007", cuts=c(-30,-2,-1,1,2,30))
plotNullres$legend$inside$args$key$points$cex<-c(1,1,1,1,1) 
print(plotNullres) #plot standardized residuals through space for one year
  #standardized residuals exhibit a large spread across space and some clustering which appears to be at a scale similar to watersheds and not related to ecological boundaries
  # later year residuals seem to be more negative on average (observed value overestimated)
Nrmse<-sqrt(mean((Nres$residuals*Nres$sd)^2)) #0.915314079

###################################
#Model A Diagnostics
###################################
A<-readRDS('/home/kate/CA/data/A')

#extract precisions for random effects, convert to variance, and extract expected value (ie mean variance)
prec<-A$marginals.hyperpar$`Precision for the Gaussian observations`
A_eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-A$marginals.hyperpar$`Precision for HUC (iid component)` 
A_v00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2) 
prec2b<-A$marginals.hyperpar$`Precision for HUC (spatial component)` 
A_s00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2b) 
prec3<-A$marginals.hyperpar$`Precision for POINTID` 
A_u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)

totMLMvarN<-eijk+u00k+u0jk #1.09124
totMLMvarA<-A_eijk+A_v00k+A_u0jk #multilevel nonstructured variance (0.543445) reduced

#compute random intercept
  #b0jk = b000 + v00k + u0jk

A_b000<-inla.rmarginal(1000,marg=A$marginals.fixed$`(Intercept)`) #random draw for the marginal distribution of the mean effect (intercept) 

#setup id sequence needed for random draw loop
id2<-as.data.frame(unique(ds$POINTID))#list of unique POINTIDs (level 2 units)
n2<- seq(1,length(id2$`unique(ds$POINTID)`))#count of level 2 units
id2<-cbind(n2,id2) 
ds<-left_join(ds,id2, by=c("POINTID"= "unique(ds$POINTID)")) #create a column in the dataset that links the original POINTID to the simplified sequence

#build a matrix of level 2 and level 3 combined random effects
A_u0jk<-matrix(NA,1000,length(n2))#an empty matrix to hold the random effects of POINTID and HUC on the intercept term b0jk
A_U<-matrix(NA,1000,length(n2)) #empty matrix to hold the level 2 and level 3 random effects
for(i in 1:length(n2)){
  A_u0jk[ ,i]<-inla.rmarginal(1000,marg=A$marginals.random$POINTID[[i]]) #random draw from the marginal distribtuion of the level 2 RE
  n3=first(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located                       
  A_v00k<-inla.rmarginal(1000,marg=A$marginals.random$HUC[1:length(HUC)][[n3]]) #random draw from the marginal distribution of the level 3 RE for the group in which the level 2 unit is nested
  A_U[ ,i]<- A_u0jk[ ,i] + A_v00k #add the appropriate level 3 random effects to the level 2 effects
}

A_beta0jk<-A_b000 + A_U #compute beta0jk
A_beta0jk_quartiles<-t(apply(A_beta0jk, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
A_beta0jk_q<-as.data.frame(A_beta0jk_quartiles)
A_beta0jk_q<-mutate(A_beta0jk_q, ID=row.names(A_beta0jk_q))
A_beta0jk_q$ID<-as.numeric(A_beta0jk_q$ID)

#plot the random intercept (beta0jk) without spatial effects
ds_dup<-ds
ds_dup<-left_join(ds, A_beta0jk_q, by=c("n2"= "ID"))
coordinates(ds_dup)= c("lon", "lat")
Aplot<-ds_dup[ ,74:77]
colnames(Aplot@data)<-c("n2","Lower", "Median", "Upper")
Aplot@proj4string <- CRS('+init=epsg:3310')
spplot(Aplot, c("Lower","Median","Upper"), main = "A: Posterior Median and 95% Credibility Interval of Random Intercept", cex =0.1)
plotA<- spplot(Aplot, c("Median"), key.space= list(x=0.5, y=0.95, corner=c(0,1)), main = "A: Posterior Median of Random Intercept", cex =0.2)
plotA$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
print(plotA)
  #shows similar spatial clustering as seen in model N, where clusters coincide largely with eco regions

#extract and plot the spatial random effect
A_s00k<-matrix(NA,1000,length(HUC))#an empty matrix to hold the spatially structured random effects of HUC on the intercept term b0jk
for(i in 1:length(HUC)){
  A_s00k[ ,i]<-inla.rmarginal(1000,marg=A$marginals.random$HUC[(length(HUC)+1):(2*length(HUC))][[i]]) #random draw from the marginal distribution of the level 3 spatially structured RE
}
A_s00k_quartiles<-t(apply(A_s00k, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
A_s00k_q<-as.data.frame(A_s00k_quartiles)
A_s00k_q<-mutate(A_s00k_q, ID=row.names(A_s00k_q))
A_s00k_q$ID<-as.numeric(A_s00k_q$ID)
 
shp_dup<-shp
shp_dup@data<-full_join(shp_dup@data, A_s00k_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
spplot(shp_dup, "Median", main = "A: Posterior Median of Watershed Spatial Effect")#plot the  median intercept
spplot(shp_dup, c("Lower","Median","Upper"), main = "A: Posterior Median and 95% Credibility Interval of Watershed Spatial Effect")
  # shows smaller regionalized clustering
  
#plot the residuals
Ares<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
Ares<-cbind(Ares,A$summary.fitted.values$mean,A$summary.fitted.values$sd)
colnames(Ares)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
Ares<-mutate(Ares, residuals = (px_tvp-mean)/sd)
plot(Ares$year,Ares$residuals, cex=0.5) #plot standardized residuals through time
  #similar to model N

Armse<-sqrt(mean((Ares$residuals*Ares$sd)^2)) #0.9152907 slightly smaller RMSE than the iid null model (0.915314079), variation at 5th decimal place

coordinates(Ares)<- ~lon+lat
Ares@proj4string <- CRS('+init=epsg:3310')
Ares@data<-Ares@data[Ares@data$year==0,]
plotAres<-spplot(Ares[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="A: Standardized \nResiduals: 2007", cuts=c(-30,-2,-1,1,2,30))
plotAres$legend$inside$args$key$points$cex<-c(1,1,1,1,1) 
print(plotAres) #plot standardized residuals through space for one year
  #generally same as model N


###################################
#Model B Diagnostics
###################################
B<-readRDS('/home/kate/CA/data/B')

#extract precisions for random effects, convert to variance, and extract expected value (ie mean variance)
prec<-B$marginals.hyperpar$`Precision for the Gaussian observations`
B_eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-B$marginals.hyperpar$`Precision for HUC (iid component)` 
B_v00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2) 
prec2b<-B$marginals.hyperpar$`Precision for HUC (spatial component)` 
B_s00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2b) #spatial effect reduced from  (0.65 for A to 0.55 for B)
prec3<-B$marginals.hyperpar$`Precision for POINTID` 
B_u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)

totMLMvarA<-A_eijk+A_v00k+A_u0jk #(0.543445)
totMLMvarB<-B_eijk+B_v00k+B_u0jk # (0.510828) multilevel nonstructured variance reduced


#compute random intercept
#b0jk = b000 + v00k + u0jk

B_b000<-inla.rmarginal(1000,marg=B$marginals.fixed$`(Intercept)`) #random draw for the marginal distribution of the mean effect (intercept) 

#setup id sequence needed for random draw loop
id2<-as.data.frame(unique(ds$POINTID))#list of unique POINTIDs (level 2 units)
n2<- seq(1,length(id2$`unique(ds$POINTID)`))#count of level 2 units
id2<-cbind(n2,id2) 
ds<-left_join(ds,id2, by=c("POINTID"= "unique(ds$POINTID)")) #create a column in the dataset that links the original POINTID to the simplified sequence


#build a matrix of level 2 and level 3 combined random effects
B_u0jk<-matrix(NA,1000,length(n2))#an empty matrix to hold the random effects of POINTID and HUC on the intercept term b0jk
B_U<-matrix(NA,1000,length(n2)) #empty matrix to hold the level 2 and level 3 random effects
for(i in 1:length(n2)){
  B_u0jk[ ,i]<-inla.rmarginal(1000,marg=B$marginals.random$POINTID[[i]]) #random draw from the marginal distribtuion of the level 2 RE
  n3=first(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located                       
  B_v00k<-inla.rmarginal(1000,marg=B$marginals.random$HUC[1:length(HUC)][[n3]]) #random draw from the marginal distribution of the level 3 RE for the group in which the level 2 unit is nested
  B_U[ ,i]<- B_u0jk[ ,i] + B_v00k #add the appropriate level 3 random effects to the level 2 effects
}

B_beta0jk<-B_b000 + B_U #compute beta0jk
B_beta0jk_quartiles<-t(apply(B_beta0jk, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
B_beta0jk_q<-as.data.frame(B_beta0jk_quartiles)
B_beta0jk_q<-mutate(B_beta0jk_q, ID=row.names(B_beta0jk_q))
B_beta0jk_q$ID<-as.numeric(B_beta0jk_q$ID)

#plot the random intercept (beta0jk) without spatial effects
ds_dup<-ds
ds_dup<-left_join(ds, B_beta0jk_q, by=c("n2"= "ID"))
coordinates(ds_dup)= c("lon", "lat")
Bplot<-ds_dup[ ,74:77]
colnames(Bplot@data)<-c("n2","Lower", "Median", "Upper")
Bplot@proj4string <- CRS('+init=epsg:3310')
spplot(Bplot, c("Lower","Median","Upper"), main = "B: Posterior Median and 95% Credibility Interval of Random Intercept", cex =0.1)
plotB<- spplot(Bplot, c("Median"), key.space= list(x=0.5, y=0.95, corner=c(0,1)), main = "B: Posterior Median of Random Intercept", cex =0.2)
plotB$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
print(plotB)
  #different appearance (maybe due to color breaks),but similar trend as model A

#extract and plot the spatial random effect
B_s00k<-matrix(NA,1000,length(HUC))#an empty matrix to hold the spatially structured random effects of HUC on the intercept term b0jk
for(i in 1:length(HUC)){
  B_s00k[ ,i]<-inla.rmarginal(1000,marg=B$marginals.random$HUC[(length(HUC)+1):(2*length(HUC))][[i]]) #random draw from the marginal distribution of the level 3 spatially structured RE
}
B_s00k_quartiles<-t(apply(B_s00k, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
B_s00k_q<-as.data.frame(B_s00k_quartiles)
B_s00k_q<-mutate(B_s00k_q, ID=row.names(B_s00k_q))
B_s00k_q$ID<-as.numeric(B_s00k_q$ID)

shp_dup<-shp
shp_dup@data<-full_join(shp_dup@data, B_s00k_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
spplot(shp_dup, "Median", main = "B: Posterior Median of Watershed Spatial Effect")#plot the  median intercept
spplot(shp_dup, c("Lower","Median","Upper"), main = "B: Posterior Median and 95% Credibility Interval of Watershed Spatial Effect")
  #similar trend as model A

#plot the residuals
Bres<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
Bres<-cbind(Bres,B$summary.fitted.values$mean,B$summary.fitted.values$sd)
colnames(Bres)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
Bres<-mutate(Bres, residuals = (px_tvp-mean)/sd)
plot(Bres$year,Bres$residuals, cex=0.5) #plot standardized residuals through time
  #similar trend as A but tightened range 

Brmse<-sqrt(mean((Bres$residuals*Bres$sd)^2)) #0.425211671 smaller RMSE than model A

coordinates(Bres)<- ~lon+lat
Bres@proj4string <- CRS('+init=epsg:3310')
Bres@data<-Bres@data[Bres@data$year==0,]
plotBres<-spplot(Bres[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="B: Standardized \nResiduals: 2007", cuts=c(-30,-2,-1,1,2,30))
plotBres$legend$inside$args$key$points$cex<-c(1,1,1,1,1) 
print(plotBres) #plot standardized residuals through space for one year
  #similar to A


###################################
#Model D Diagnostics
###################################
D<-readRDS('/home/kate/CA/data/D')

#extract precisions for random effects, convert to variance, and extract expected value (ie mean variance)
prec<-D$marginals.hyperpar$`Precision for the Gaussian observations`
D_eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-D$marginals.hyperpar$`Precision for HUC (iid component)` 
D_v00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2) 
prec2b<-D$marginals.hyperpar$`Precision for HUC (spatial component)` 
D_s00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2b) #spatial effect (0.551) slightly increased over model B
prec3<-D$marginals.hyperpar$`Precision for POINTID` 
D_u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)

totMLMvarA<-A_eijk+A_v00k+A_u0jk 
totMLMvarB<-B_eijk+B_v00k+B_u0jk #multilevel nonstructured variance reduced
totMLMvarD<-D_eijk+D_v00k+D_u0jk #multilevel nonstructured variance (0.5177735) reduced a tiny bit

#compute random intercept
#b0jk = b000 + v00k + u0jk

D_b000<-inla.rmarginal(1000,marg=D$marginals.fixed$`(Intercept)`) #random draw for the marginal distribution of the mean effect (intercept) 

# #setup id sequence needed for random draw loop
# id2<-as.data.frame(unique(ds$POINTID))#list of unique POINTIDs (level 2 units)
# n2<- seq(1,length(id2$`unique(ds$POINTID)`))#count of level 2 units
# id2<-cbind(n2,id2) 
# ds<-left_join(ds,id2, by=c("POINTID"= "unique(ds$POINTID)")) #create a column in the dataset that links the original POINTID to the simplified sequence


#build a matrix of level 2 and level 3 combined random effects
D_u0jk<-matrix(NA,1000,length(n2))#an empty matrix to hold the random effects of POINTID and HUC on the intercept term b0jk
D_U<-matrix(NA,1000,length(n2)) #empty matrix to hold the level 2 and level 3 random effects
for(i in 1:length(n2)){
  D_u0jk[ ,i]<-inla.rmarginal(1000,marg=D$marginals.random$POINTID[[i]]) #random draw from the marginal distribtuion of the level 2 RE
  n3=first(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located                       
  D_v00k<-inla.rmarginal(1000,marg=D$marginals.random$HUC[1:length(HUC)][[n3]]) #random draw from the marginal distribution of the level 3 RE for the group in which the level 2 unit is nested
  D_U[ ,i]<- D_u0jk[ ,i] + D_v00k #add the appropriate level 3 random effects to the level 2 effects
}

D_beta0jk<-D_b000 + D_U #compute beta0jk
D_beta0jk_quartiles<-t(apply(D_beta0jk, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
D_beta0jk_q<-as.data.frame(D_beta0jk_quartiles)
D_beta0jk_q<-mutate(D_beta0jk_q, ID=row.names(D_beta0jk_q))
D_beta0jk_q$ID<-as.numeric(D_beta0jk_q$ID)

#plot the random intercept (beta0jk) without spatial effects
ds_dup<-ds
ds_dup<-left_join(ds, D_beta0jk_q, by=c("n2"= "ID"))
coordinates(ds_dup)= c("lon", "lat")
Dplot<-ds_dup[ ,74:77]
colnames(Dplot@data)<-c("n2","Lower", "Median", "Upper")
Dplot@proj4string <- CRS('+init=epsg:3310')
spplot(Dplot, c("Lower","Median","Upper"), main = "D: Posterior Median and 95% Credibility Interval of Random Intercept", cex =0.1)
plotD<- spplot(Dplot, c("Median"), key.space= list(x=0.5, y=0.95, corner=c(0,1)), main = "D: Posterior Median of Random Intercept", cex =0.2)
plotD$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
print(plotD)
#note

#extract and plot the spatial random effect
D_s00k<-matrix(NA,1000,length(HUC))#an empty matrix to hold the spatially structured random effects of HUC on the intercept term b0jk
for(i in 1:length(HUC)){
  D_s00k[ ,i]<-inla.rmarginal(1000,marg=D$marginals.random$HUC[(length(HUC)+1):(2*length(HUC))][[i]]) #random draw from the marginal distribution of the level 3 spatially structured RE
}
D_s00k_quartiles<-t(apply(D_s00k, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
D_s00k_q<-as.data.frame(D_s00k_quartiles)
D_s00k_q<-mutate(D_s00k_q, ID=row.names(D_s00k_q))
D_s00k_q$ID<-as.numeric(D_s00k_q$ID)

shp_dup<-shp
shp_dup@data<-full_join(shp_dup@data, D_s00k_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
spplot(shp_dup, "Median", main = "D: Posterior Median of Watershed Spatial Effect")#plot the  median intercept
spplot(shp_dup, c("Lower","Median","Upper"), main = "D: Posterior Median and 95% Credibility Interval of Watershed Spatial Effect")
#note

#plot the residuals
Dres<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
Dres<-cbind(Dres,D$summary.fitted.values$mean,D$summary.fitted.values$sd)
colnames(Dres)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
Dres<-mutate(Dres, residuals = (px_tvp-mean)/sd)
plot(Dres$year,Dres$residuals, cex=0.5) #plot standardized residuals through time
#pertty much the same as B

Drmse<-sqrt(mean((Dres$residuals*Dres$sd)^2)) #0.42515577 smaller RMSE than model B

coordinates(Dres)<- ~lon+lat
Dres@proj4string <- CRS('+init=epsg:3310')
Dres@data<-Dres@data[Dres@data$year==7,]
plotDres<-spplot(Dres[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="D: Standardized \nResiduals: 2014", cuts=c(-30,-2,-1,1,2,30))
plotDres$legend$inside$args$key$points$cex<-c(1,1,1,1,1) 
print(plotDres) #plot standardized residuals through space for one year
#note

###################################
#Model E Diagnostics
###################################
E<-readRDS('/home/kate/CA/data/E')

#extract precisions for random effects, convert to variance, and extract expected value (ie mean variance)
prec<-E$marginals.hyperpar$`Precision for the Gaussian observations`
E_eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-E$marginals.hyperpar$`Precision for HUC (iid component)` 
E_v00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2) 
prec2b<-E$marginals.hyperpar$`Precision for HUC (spatial component)` 
E_s00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2b) #spatial effect (0.565) increased over D, but still smaller than A
prec3<-E$marginals.hyperpar$`Precision for POINTID` 
E_u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)


totMLMvarE<-E_eijk+E_v00k+E_u0jk #multilevel nonstructured variance reduced from model D (0.510777 to 0.462350671)

#compute random intercept
#b0jk = b000 + v00k + u0jk

E_b000<-inla.rmarginal(1000,marg=E$marginals.fixed$`(Intercept)`) #random draw for the marginal distribution of the mean effect (intercept) 

# #setup id sequence needed for random draw loop
# id2<-as.data.frame(unique(ds$POINTID))#list of unique POINTIDs (level 2 units)
# n2<- seq(1,length(id2$`unique(ds$POINTID)`))#count of level 2 units
# id2<-cbind(n2,id2) 
# ds<-left_join(ds,id2, by=c("POINTID"= "unique(ds$POINTID)")) #create a column in the dataset that links the original POINTID to the simplified sequence


#build a matrix of level 2 and level 3 combined random effects
E_u0jk<-matrix(NA,1000,length(n2))#an empty matrix to hold the random effects of POINTID and HUC on the intercept term b0jk
E_U<-matrix(NA,1000,length(n2)) #empty matrix to hold the level 2 and level 3 random effects
for(i in 1:length(n2)){
  E_u0jk[ ,i]<-inla.rmarginal(1000,marg=E$marginals.random$POINTID[[i]]) #random draw from the marginal distribtuion of the level 2 RE
  n3=first(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located                       
  E_v00k<-inla.rmarginal(1000,marg=E$marginals.random$HUC[1:length(HUC)][[n3]]) #random draw from the marginal distribution of the level 3 RE for the group in which the level 2 unit is nested
  E_U[ ,i]<- E_u0jk[ ,i] + E_v00k #add the appropriate level 3 random effects to the level 2 effects
}

E_beta0jk<-E_b000 + E_U #compute beta0jk
E_beta0jk_quartiles<-t(apply(E_beta0jk, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
E_beta0jk_q<-as.data.frame(E_beta0jk_quartiles)
E_beta0jk_q<-mutate(E_beta0jk_q, ID=row.names(E_beta0jk_q))
E_beta0jk_q$ID<-as.numeric(E_beta0jk_q$ID)

#plot the random intercept (beta0jk) without spatial effects
ds_dup<-ds
ds_dup<-left_join(ds, E_beta0jk_q, by=c("n2"= "ID"))
coordinates(ds_dup)= c("lon", "lat")
Eplot<-ds_dup[ ,74:77]
colnames(Eplot@data)<-c("n2","Lower", "Median", "Upper")
Eplot@proj4string <- CRS('+init=epsg:3310')
spplot(Eplot, c("Lower","Median","Upper"), main = "E: Posterior Median and 95% Credibility Interval of Random Intercept", cex =0.1)
plotE<- spplot(Eplot, c("Median"), key.space= list(x=0.5, y=0.95, corner=c(0,1)), main = "E: Posterior Median of Random Intercept", cex =0.2)
plotE$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
print(plotE)
#note

#extract and plot the spatial random effect
E_s00k<-matrix(NA,1000,length(HUC))#an empty matrix to hold the spatially structured random effects of HUC on the intercept term b0jk
for(i in 1:length(HUC)){
  E_s00k[ ,i]<-inla.rmarginal(1000,marg=E$marginals.random$HUC[(length(HUC)+1):(2*length(HUC))][[i]]) #random draw from the marginal distribution of the level 3 spatially structured RE
}
E_s00k_quartiles<-t(apply(E_s00k, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975)))) #compute quartiles for b0jk
E_s00k_q<-as.data.frame(E_s00k_quartiles)
E_s00k_q<-mutate(E_s00k_q, ID=row.names(E_s00k_q))
E_s00k_q$ID<-as.numeric(E_s00k_q$ID)

shp_dup<-shp
shp_dup@data<-full_join(shp_dup@data, E_s00k_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
spplot(shp_dup, "Median", main = "E: Posterior Median of Watershed Spatial Effect")#plot the  median intercept
spplot(shp_dup, c("Lower","Median","Upper"), main = "E: Posterior Median and 95% Credibility Interval of Watershed Spatial Effect")
#note

#plot the residuals
Eres<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
Eres<-cbind(Eres,E$summary.fitted.values$mean,E$summary.fitted.values$sd)
colnames(Eres)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
Eres<-mutate(Eres, residuals = (px_tvp-mean)/sd)
plot(Eres$year,Eres$residuals, cex=0.5) #plot standardized residuals through time
#actually slightly bigger spread within years, with trend remaining about the same

Ermse<-sqrt(mean((Eres$residuals*Eres$sd)^2)) #0.368585 smaller RMSE than model D

coordinates(Eres)<- ~lon+lat
Eres@proj4string <- CRS('+init=epsg:3310')
Eres@data<-Eres@data[Eres@data$year==0,]
plotEres<-spplot(Eres[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="E: Standardized \nResiduals: 2007", cuts=c(-30,-2,-1,1,2,30))
plotEres$legend$inside$args$key$points$cex<-c(1,1,1,1,1) 
print(plotEres) #plot standardized residuals through space for one year
#note


################################
#################################
###ENTERING OLD SCRIPT TERRITORY###
###############################
######################################

#Area Map
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
#Raw data plots
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


###These are pretty much flat - negligible time variance in watershed level water rights vars###

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
b0<-stlm3e$summary.fixed[1,] #overall intercept
b1<-stlm3e$summary.fixed[12,] #main effect of the focal predictor, rip in this case
b2<-stlm3e$summary.fixed[2,]#main effect of the moderator, year
b3<-stlm3e$summary.fixed[15,] #interaction effect, rip:year

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
b1s<-inla.rmarginal(1000,marg=stlm3e$marginals.fixed$rip_perc)
b3s<-inla.rmarginal(1000, marg=stlm3e$marginals.fixed$`year:rip_perc`)
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
b0<-stlm3e$summary.fixed[1,] #overall intercept
b1<-stlm3e$summary.fixed[13,] #main effect of the focal predictor, rip in this case
b2<-stlm3e$summary.fixed[2,]#main effect of the moderator, year
b3<-stlm3e$summary.fixed[16,] #interaction effect, rip:year

plot(focal, ((b0[,1] + b2[,1]*moderator[2])+(b1[,1] + b3[,1]*moderator[2])*focal), type="l", ylim=c(0.42,0.51), xlab="Percent Appropriative Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Appropriative Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1[,1] + b3[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1[,1] + b3[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[8])+(b1[,1] + b3[,1]*moderator[8])*focal), col="red")
legend(0.7,0.51, c("1", "3", "5", "7"),title ="Years of Drought", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","green","blue", "red"))

mod<-as.matrix(moderator[c(1,2,4,6,8)])
b1s<-inla.rmarginal(1000,marg=stlm3e$marginals.fixed$approp_perc)
b3s<-inla.rmarginal(1000, marg=stlm3e$marginals.fixed$`year:approp_perc`)
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
b0<-stlm3e$summary.fixed[1,] #overall intercept
b1<-stlm3e$summary.fixed[14,] #main effect of the focal predictor, rip in this case
b2<-stlm3e$summary.fixed[2,]#main effect of the moderator, year
b3<-stlm3e$summary.fixed[17,] #interaction effect, rip:year

plot(focal, ((b0[,1] + b2[,1]*moderator[2])+(b1[,1] + b3[,1]*moderator[2])*focal), type="l", ylim=c(0.42,0.54), xlab="Percent Pre-1914 Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Pre-1914 Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1[,1] + b3[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1[,1] + b3[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[8])+(b1[,1] + b3[,1]*moderator[8])*focal), col="red")
legend(0.7,0.54, c("1", "3", "5", "7"),title ="Years of Drought", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("black","green","blue", "red"))

mod<-as.matrix(moderator[c(1,2,4,6,8)])
b1s<-inla.rmarginal(1000,marg=stlm3e$marginals.fixed$pre_perc)
b3s<-inla.rmarginal(1000, marg=stlm3e$marginals.fixed$`year:pre_perc`)
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
b0p<-stlm3e$summary.fixed[1,] #overall intercept
b1p<-stlm3e$summary.fixed[14,] #main effect of the focal predictor, pre in this case
b2p<-stlm3e$summary.fixed[2,]#main effect of the moderator, year
b3p<-stlm3e$summary.fixed[17,] #interaction effect, pre:year
b0a<-stlm3e$summary.fixed[1,] #overall intercept
b1a<-stlm3e$summary.fixed[13,] #main effect of the focal predictor, approp in this case
b2a<-stlm3e$summary.fixed[2,]#main effect of the moderator, year
b3a<-stlm3e$summary.fixed[16,] #interaction effect, approp:year
b0r<-stlm3e$summary.fixed[1,] #overall intercept
b1r<-stlm3e$summary.fixed[12,] #main effect of the focal predictor, rip in this case
b2r<-stlm3e$summary.fixed[2,]#main effect of the moderator, year
b3r<-stlm3e$summary.fixed[15,] #interaction effect, rip:year
r1<-stlm3c$summary.random$year1[2,1] #random time effect year 1 of 0-7
r4<-stlm3c$summary.random$year1[5,1]
r7<-stlm3c$summary.random$year1[8,1]
r<-stlm3c$summary.random$year1[ ,1]

plot(focal, ((b0p[,1] + b2p[,1]*moderator[1] )+(b1p[,1] + b3p[,1]*moderator[1])*focal), type="l", col="red", lwd=2,ylim=c(-0.5,0.5), xlab="Percent Pre-1914 Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(focal, ((b0p[,1] + b2p[,1]*moderator[4] )+(b1p[,1] + b3p[,1]*moderator[4])*focal), col="red", lty=2, lwd=2)
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
plot(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[2]), type="l", col="red",ylim=c(-.5,-0.1), xlab="Year of Drought", 
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

#This plot flips the x axis and includes the year fixed effect
plot(moderator, ((b0p[,1] + b2p[,1]*moderator + r[moderator+1])+(b1p[,1] + b3p[,1]*moderator)*focal[2]), type="l", col="red",ylim=c(-.5,-0.1), xlab="Year of Drought", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Water Rights as a Function of Time",cex.main = 1 ) #conditional regression on tvp as fcn of time
lines(moderator, ((b0p[,1] + b2p[,1]*moderator + r[moderator+1])+(b1p[,1] + b3p[,1]*moderator)*focal[3]), col="red", lty=2)
lines(moderator, ((b0p[,1] + b2p[,1]*moderator + r[moderator+1])+(b1p[,1] + b3p[,1]*moderator)*focal[4]), col="red", lty=4)
lines(moderator, ((b0a[,1] + b2a[,1]*moderator + r[moderator+1])+(b1a[,1] + b3a[,1]*moderator)*focal[2]), col="black")
lines(moderator, ((b0a[,1] + b2a[,1]*moderator + r[moderator+1])+(b1a[,1] + b3a[,1]*moderator)*focal[3]), col="black", lty=2)
lines(moderator, ((b0a[,1] + b2a[,1]*moderator + r[moderator+1])+(b1a[,1] + b3a[,1]*moderator)*focal[4]), col="black", lty=4)
lines(moderator, ((b0r[,1] + b2r[,1]*moderator + r[moderator+1])+(b1r[,1] + b3r[,1]*moderator)*focal[2]), col="blue")
lines(moderator, ((b0r[,1] + b2r[,1]*moderator + r[moderator+1])+(b1r[,1] + b3r[,1]*moderator)*focal[3]), col="blue", lty=2)
lines(moderator, ((b0r[,1] + b2r[,1]*moderator + r[moderator+1])+(b1r[,1] + b3r[,1]*moderator)*focal[4]), col="blue", lty=4)
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
res<-cbind(res,p1fv$summary.fitted.values$mean,stlm3cfv$summary.fitted.values$sd)
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



###################################################
#Random wr effect plot 
#########################################

#for model that is same as stlm3c but with water rights random varying across eyars according to rw1
plot(stlm4b$summary.random$year5[ ,1],stlm4b$summary.random$year5[ ,2], ylim=c(-0.08, 0.08) ) #pre-1914 mean
points(stlm4b$summary.random$year4[ ,1],stlm4b$summary.random$year4[ ,2], col="red") #approp mean
points(stlm4b$summary.random$year3[ ,1],stlm4b$summary.random$year3[ ,2], col="blue") #rip mean
 points(stlm4b$summary.random$year5[ ,1],stlm4b$summary.random$year5[ ,4], cex=.5) #pre 1914 2.5%
 points(stlm4b$summary.random$year5[ ,1],stlm4b$summary.random$year5[ ,6], cex=.5) #pre 1914 97.5 %
points(stlm4b$summary.random$year4[ ,1],stlm4b$summary.random$year4[ ,6], cex=.5, col="red")
points(stlm4b$summary.random$year4[ ,1],stlm4b$summary.random$year4[ ,4], cex=.5, col="red")
 points(stlm4b$summary.random$year3[ ,1],stlm4b$summary.random$year3[ ,6], cex=.5, col="blue")
 points(stlm4b$summary.random$year3[ ,1],stlm4b$summary.random$year3[ ,4], cex=.5, col="blue")
 
#model that includes space effect, time effect, and space time interactions as well as randomly varying water rights (rw1 across eyars)
 plot(nonlin3$summary.random$year6[ ,1],nonlin3$summary.random$year6[ ,2], ylim=c(-0.1, 0.15), col="red" )
 points(nonlin3$summary.random$year5[ ,1],nonlin3$summary.random$year5[ ,2])
 points(nonlin3$summary.random$year4[ ,1],nonlin3$summary.random$year4[ ,2], col="blue")
 points(nonlin3$summary.random$year6[ ,1],nonlin3$summary.random$year6[ ,4], cex=.5, , col="red")
 points(nonlin3$summary.random$year6[ ,1],nonlin3$summary.random$year6[ ,6], cex=.5, , col="red")
 points(nonlin3$summary.random$year5[ ,1],nonlin3$summary.random$year5[ ,6], cex=.5)
 points(nonlin3$summary.random$year5[ ,1],nonlin3$summary.random$year5[ ,4], cex=.5)
 points(nonlin3$summary.random$year4[ ,1],nonlin3$summary.random$year4[ ,6], cex=.5, col="blue")
 points(nonlin3$summary.random$year4[ ,1],nonlin3$summary.random$year4[ ,4], cex=.5, col="blue")
 
 #model that includes space effect, time effect, and space time interactions as well as randomly varying water rights (rw1 across eyars)
 plot(fe2$summary.random$year6[ ,1],fe2$summary.random$year6[ ,2], ylim=c(-0.1, 0.15), col="red" )
 points(fe2$summary.random$year5[ ,1],fe2$summary.random$year5[ ,2])
 points(fe2$summary.random$year4[ ,1],fe2$summary.random$year4[ ,2], col="blue")
 points(fe2$summary.random$year6[ ,1],fe2$summary.random$year6[ ,4], cex=.5, , col="red")
 points(fe2$summary.random$year6[ ,1],fe2$summary.random$year6[ ,6], cex=.5, , col="red")
 points(fe2$summary.random$year5[ ,1],fe2$summary.random$year5[ ,6], cex=.5)
 points(fe2$summary.random$year5[ ,1],fe2$summary.random$year5[ ,4], cex=.5)
 points(fe2$summary.random$year4[ ,1],fe2$summary.random$year4[ ,6], cex=.5, col="blue")
 points(fe2$summary.random$year4[ ,1],fe2$summary.random$year4[ ,4], cex=.5, col="blue")