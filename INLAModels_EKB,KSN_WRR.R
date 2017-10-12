library(memisc)
library(dplyr)
library (gstat)
library(sp)
library(nlme)
library(rgdal)
library(INLA)
library(spdep)
library(raster)
library(lattice)
library(ggplot2)

setwd('C:/')

#This script reads in the full compiled California water rights dataset, builds additionallagged and binary variables, 
#builds spatial relationship objects using R-INLA, runs models in R-INLA, and prepares plots and figures for model result 
#interpretation. This script was used in the preparation of the manuscript Nelson, K. S., & Burchfield, E. K. Effects of 
#the Structure of Water Rights on Agricultural Production During Drought: A Spatiotemporal Analysis of California's 
#Central Valley. Water Resources Research.DOI: 10.1002/2017WR020666.

##############################################
#Read in datasets to use in INLA models
##############################################

###  full dataset with lat lon columns and updated wr info --> check for cleanliness 
d<- readRDS('/home/kate/CA/data/cv.dat06072017.rds') 
#plot(d)

###add water contract boundary files and reconcile projection
cvp<-shapefile('/home/kate/CA/data/CVPcontractors') #central valley project water contractor district boundaries
swp<-shapefile('/home/kate/CA/data/wdst24') #state water project water contractor boundaries
pwc<-shapefile('/home/kate/CA/data/wdpr24') #private water contractor boundaries
swp@proj4string<- CRS(proj4string(pwc)) #supply projection info for state water project shapefile (metadata indicates same proj as pwc)

cvp <- spTransform(cvp, CRS(proj4string(d)))
swp <- spTransform(swp, CRS(proj4string(d)))
pwc <- spTransform(pwc, CRS(proj4string(d)))

# plot(cvp, col="red")
# plot(swp, add=TRUE, col="blue")
# plot(pwc, add=TRUE, col="green")

#####################################
#Build water contract boundary binary
#######################################
dsub<-d
dsub<-dsub[dsub@data$year==14, ]
dsub@data<-dsub@data[ ,1:5]
loc<-over(dsub,cvp[ ,"WDNAME"])
dnew<-dsub
dnew@data <- cbind(dsub@data, loc)
names(dnew@data)[names(dnew@data) == 'WDNAME'] <- 'WD_CVP'
loc2<-over(dsub, swp[,"WDNAME"])
dnew@data<-cbind(dnew@data, loc2)
names(dnew@data)[names(dnew@data) == 'WDNAME'] <- 'WD_SWP'
loc3<-over(dsub, pwc[,"WDNAME"])
dnew@data<-cbind(dnew@data, loc3)
names(dnew@data)[names(dnew@data) == 'WDNAME'] <- 'WD_PWC'

# plot(dnew[!is.na(dnew@data$WD_CVP),])
# points(dnew[!is.na(dnew@data$WD_SWP),], col="red")
# points(dnew[!is.na(dnew@data$WD_PWC),], col="green")
# plot(cv, add=TRUE)

d_update<-d
d_update<-merge(d_update, dnew@data, by=c("POINTID"), duplicateGeom=TRUE) #join to original space-time dataset
names(d_update@data)[names(d_update@data) == 'FID.x'] <- 'FID'
names(d_update@data)[names(d_update@data) == 'GRID_CODE.x'] <- 'GRID_CODE'
names(d_update@data)[names(d_update@data) == 'huc12_id.x'] <- 'huc12_id'
names(d_update@data)[names(d_update@data) == 'AreaAcres.x'] <- 'AreaAcres'

#create binary
d_update$govwc<-is.na(d_update$WD_CVP) & is.na(d_update@data$WD_SWP)
d_update$govwc<-((1*d_update$govwc)-1)*-1  #1 for in a gov wd and 0 for not in a government water district
d_update$pwc<-!is.na(d_update$WD_PWC) 
d_update$pwc<-1*d_update$pwc #1 for in a private wd and 0 for not in 
d_update$wd<-is.na(d_update$WD_CVP) & is.na(d_update@data$WD_SWP) & is.na(d_update$WD_PWC) 
d_update$wd<-((1*d_update$wd)-1)*-1 #1 for in any water district

d<-d_update #update the original dataset

###################################################################################################
#Offset the farmland crop diversity measure,SPI, and copr type by one year and add lat/lon values
####################################################################################################

yrset<-c(7:14)
lagset<-d[d@data$year %in% yrset, c(1,11,13,15, 17)]
lagset$year<-lagset$year + 1 #want to associate a past year (2007) with the current year (2008)
d2<-d #work with a duplicate
d2<-merge(d2, lagset, by=c("POINTID"="POINTID", "year"="year"))
names(d2@data)[names(d2@data) == 'spi.x'] <- 'spi' #leave unlagged variables with original name, lagged variables will be name.y
names(d2@data)[names(d2@data) == 'fl_cd.x'] <- 'fl_cd'
names(d2@data)[names(d2@data) == 'aglulc.x'] <- 'aglulc'
d<-d2

sub<-d
sub<-sub[sub@data$year==7,] #pull out a single year or coordinates will get jumbled
sub@data$lat<-sub@coords[ ,2]
sub@data$lon<-sub@coords[ ,1]
coordset<-sub@data[ ,c(1,63,64)]
d2<-merge(d2, coordset, by =c("POINTID"="POINTID"))
d<-d2


##############################################
#Select farmland subset for analysis
##############################################

#full central valley watershed as spatial polygons data frame
#cv <- shapefile('/data/emily/WF/kate/shp/cv_huc.shp')
cv <- shapefile('/home/kate/CA/data/cv_huc.shp')
cv <- spTransform(cv, CRS(proj4string(d))) #set to same proj
proj4string(cv)==proj4string(d) #double check proj

#plot together to check validity of spatial overlay
# plot(d, cex=0.05, col="blue")
# lines(cv)

#remove uneeded info from cv
cv<-cv[ ,8]

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
#order by huc then time like they do in the Ohio example chapter 7 of r-inla book --> works for bym model
ds <- ds[order(ds$HUC,ds$year),]

#create neighbors list from polygon object
temp <- poly2nb(shp, queen=FALSE) #neighbors must share more than one point at boundary
#temp <- poly2nb(shp, queen=T)

#convert to a sparse matrix to reduce memory
H.adj <- nb2mat(temp, style ="B") #convert neighbor list to binary coded neighbor weights matrix
H.adj <-as(H.adj, "dgTMatrix") #convert to a sparse sytle matrix

#plot spatial object
image(inla.graph2matrix(H.adj), xlab="", ylab="")

########################################
#Prep data for analysis
#######################################
# dataset<-ds[,c(1,5,7:8,10:15,17:36, 38:50, 57:67)] #dataset for Supporting Information
# write.csv(dataset,'dataset_SI_WRR.csv')

y=ds$px_tvp #dependent variable
ds$year <- ds$year - 7 #center time at 2007 = 0
ds<-ds[ds$year>2,] #remove years where digitized reports for Statement of Water Use and Diversion are not consistently available (prior to 2009/2010)

#build duplicate indexes for model building
ds <- mutate(ds, year1=year, year2=year, year3=year) #add duplicate year grouping index
ds<-mutate(ds,POINTID1=POINTID)#add duplicate pixel grouping index
ds <- mutate(ds, HUC1=HUC,HUC2=HUC, HUC3=HUC) #add duplicates of HUC index for random slope grouping
ds$agluid<-as.numeric(ds$aglulc)#landuse id var for function fields

#save a copy of data pre-standardization for later use
dsdup<-ds
ds[mapply(is.infinite, ds)] <- NA #get rid of Inf

#add nonlinear spi terms
ds$spi2<-ds$spi*ds$spi
ds$spi3<-ds$spi*ds$spi*ds$spi
ds$spi4<-ds$spi*ds$spi*ds$spi*ds$spi

#standardize vars to ease interpretation
ds$gw<- (ds$gw - mean(ds$gw, na.rm=T))/sd(ds$gw, na.rm=T)
ds$spi<- (ds$spi - mean(ds$spi, na.rm=T))/sd(ds$spi, na.rm=T) 
ds$spi.y<- (ds$spi.y - mean(ds$spi.y, na.rm=T))/sd(ds$spi.y, na.rm=T) 
ds$fl_cd.y<-(ds$fl_cd.y - mean(ds$fl_cd.y, na.rm=TRUE))/sd(ds$fl_cd.y, na.rm=TRUE)

ds$wr_dens<- (ds$wr_dens - mean(ds$wr_dens, na.rm=T))/sd(ds$wr_dens,na.rm=T)                             
ds$pre_perc<- (ds$pre_perc - mean(ds$pre_perc, na.rm=T))/sd(ds$pre_perc, na.rm=T)
ds$rip_perc<- (ds$rip_perc - mean(ds$rip_perc, na.rm=T))/sd(ds$rip_perc, na.rm=T)
ds$approp_perc<- (ds$approp_perc - mean(ds$approp_perc, na.rm=T))/sd(ds$approp_perc, na.rm=T)

ds$ag_dens<- (ds$ag_dens - mean(ds$ag_dens, na.rm=T))/sd(ds$ag_dens,na.rm=T) 
ds$agpre_perc<- (ds$agpre_perc - mean(ds$agpre_perc, na.rm=T))/sd(ds$agpre_perc, na.rm=T)
ds$agrip_perc<- (ds$agrip_perc - mean(ds$agrip_perc, na.rm=T))/sd(ds$agrip_perc, na.rm=T)
ds$agapprop_perc<- (ds$agapprop_perc - mean(ds$agapprop_perc, na.rm=T))/sd(ds$agapprop_perc, na.rm=T)

ds$ag_perc<-(ds$ag_perc - mean(ds$ag_perc, na.rm=T))/sd(ds$ag_perc,na.rm=T) 
ds$dom_perc<- (ds$dom_perc - mean(ds$dom_perc, na.rm=T))/sd(ds$dom_perc,na.rm=T) 
ds$ind_perc<- (ds$ind_perc - mean(ds$ind_perc, na.rm=T))/sd(ds$ind_perc,na.rm=T) 
ds$fish_perc<- (ds$fish_perc - mean(ds$fish_perc, na.rm=T))/sd(ds$fish_perc,na.rm=T) 
ds$rec_perc<- (ds$rec_perc - mean(ds$rec_perc, na.rm=T))/sd(ds$rec_perc,na.rm=T) 

ds$anyagwr<-ds$ag_cnt
ds$anyagwr[ds$anyagwr != 0] <- 1 #watershed without ag water rights  and all other 1

ds$spi2<- (ds$spi2 - mean(ds$spi2, na.rm=T))/sd(ds$spi2, na.rm=T) 
ds$spi3<- (ds$spi3 - mean(ds$spi3, na.rm=T))/sd(ds$spi3, na.rm=T) 
ds$spi4<- (ds$spi4 - mean(ds$spi4, na.rm=T))/sd(ds$spi4, na.rm=T)

ds$px_tvp<- (ds$px_tvp - mean(ds$px_tvp, na.rm=T))/sd(ds$px_tvp, na.rm=T)
y<-ds$px_tvp


#Dependent var for Bernoulli likelihood, logistic model 
ds_B<-ds[ds$year>0, ]
y_B<-ds_B$aglulc
y_B[y_B!=1 ]<-0
y_B<-as.character(y_B)
y_B<-as.data.frame(as.numeric(y_B)) #barren and fallow is 1, all other landuse cats are 0


##############
### Models ###
##############

##Null 3 level ## 
#yijk = b000 + u00k + u0jk + eijk 
NULL_MOD <-y ~ 1 + f(HUC, model="iid") + f(POINTID, model = "iid")
OUT <- inla (NULL_MOD, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/N2')
#N<-readRDS('/home/kate/CA/data/final/N')


##Space-time model ## 
#yijk = b000 + u00k + u0jk + ui00 + eijk, where u00k = v00k + s00k (area specific and spatially structured randomness)
SpT <-y ~ 1 + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)  + factor(year)
OUT <- inla (SpT, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/SpT')
SpT<-readRDS('/home/kate/CA/data/final/SpT')


##Model without spatial### 
model7c <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi + 
  gw + ag_dens + ag_perc + factor(aglulc) + factor(year)
OUT <- inla (model7c, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/model7c2')
#model7c<-readRDS('/home/kate/CA/data/final/model7c')
rm(OUT)

##Model without spatial###  remove ag_perc control
model7d <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi + 
  gw + ag_dens  + factor(aglulc) + factor(year)
OUT <- inla (model7d, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/model7d2')
#model7d<-readRDS('/home/kate/CA/data/final/model7d')
rm(OUT)

##Model without spatial###  remove gw control
model7e <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi + 
  ag_dens + ag_perc  + factor(aglulc) + factor(year)
OUT <- inla (model7e, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/model7e2')
#model7e<-readRDS('/home/kate/CA/data/final/model7e2')
rm(OUT)

##Model without spatial###  remove ag_dens control
model7f <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi + 
  gw + ag_perc  + factor(aglulc) + factor(year)
OUT <- inla (model7f, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/model7f2')
#model7f<-readRDS('/home/kate/CA/data/final/model7f')
rm(OUT)

##Model full ###
modelfulla <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi +
  gw + ag_dens + ag_perc + factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (modelfulla, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelfulla2')
#modelfulla2<-readRDS('/home/kate/CA/data/final/modelfulla2')

##Model full ### remove ag_perc
modelfullb <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi +
  gw + ag_dens + factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (modelfullb, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelfullb2')
#modelfullb<-readRDS('/home/kate/CA/data/final/modelfullb')

##Model full ### remove ag_dens
modelfullc <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi +
  gw  + ag_perc + factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (modelfullc, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelfullc2')
#modelfullc<-readRDS('/home/kate/CA/data/final/modelfullc')

###Model full ### remove gw
modelfulld <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi +
  ag_perc + ag_dens + factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (modelfulld, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelfulld2b')
#modelfulld<-readRDS('/home/kate/CA/data/final/modelfulld2')

##Model full ### without riparian
modelfulle <-y ~ 1 +  agpre_perc+ agapprop_perc + spi + agpre_perc*spi + agapprop_perc*spi +
  gw + ag_dens + ag_perc + factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (modelfulle, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelfulle')
#modelfulle<-readRDS('/home/kate/CA/data/final/modelfulle')

#Model full ### add contract districts
modelfullf <-y ~ 1 +  agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi +
  gw + ag_dens + ag_perc + factor(wd)+ factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (modelfullf, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelfullfb')
#modelfulle<-readRDS('/home/kate/CA/data/final/modelfulle')

##Model FE### Use huc fixed effects
modelFE <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi + 
  gw + ag_dens  + ag_perc + factor(aglulc) + factor(year) + factor(HUC)
OUT <- inla (modelFE, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/modelFEfullb')
modelFE<-readRDS('/home/kate/CA/data/final/modelFEfull')

##quadratic
MODq <-y ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi + spi2 + agrip_perc*spi + agpre_perc*spi + agapprop_perc*spi + 
  gw + ag_dens  + ag_perc + factor(aglulc) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
OUT <- inla (MODq, family = "gaussian", data=ds, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(OUT, file='/home/kate/CA/data/final/quad')
quad<-readRDS('/home/kate/CA/data/final/quad')
# summary(quad)
# s<-seq(-3,3)
# ys<-s*quad$summary.fixed[7,4]+s^2*quad$summary.fixed[8,4]
# plot(s,ys)

##Bernoulli model##
#Bernoulli likelihood of landuse category being barren&fallow

#spacetime
Bspt <-y_B ~ 1 + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (Bspt, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bspt')
#B1b<-readRDS('/home/kate/CA/data/final/B1b')
rm(B_OUT)

#no spatial
B6b <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y +
  gw + ag_dens + ag_perc + fl_cd.y + factor(year)
B_OUT <- inla (B6b, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/B6b2b')
#B6b<-readRDS('/home/kate/CA/data/final/B6b2')
rm(B_OUT)

#no spatial, remove ag_dens
B6c <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y +
  gw  + ag_perc + fl_cd.y + factor(year)
B_OUT <- inla (B6c, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/B6c2b')
#B6b<-readRDS('/home/kate/CA/data/final/B6b')
rm(B_OUT)

#no spatial, remove ag_perc
B6d <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y +
  gw  + ag_dens + fl_cd.y + factor(year)
B_OUT <- inla (B6d, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/B6d2b')
#B6b<-readRDS('/home/kate/CA/data/final/B6b')
rm(B_OUT)

#no spatial, remove fl.cd
B6e<-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y +
  gw  + ag_dens + ag_perc + factor(year)
B_OUT <- inla (B6e, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/B6e2b')
#B6b<-readRDS('/home/kate/CA/data/final/B6b')
rm(B_OUT)

#no spatial, remove gw
B6f <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y +
  ag_dens + ag_perc + fl_cd.y + factor(year)
B_OUT <- inla (B6f, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/B6f2b')
#B6f2<-readRDS('/home/kate/CA/data/final/B6f2')
rm(B_OUT)

#full model
B_MODb <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  gw + ag_dens + ag_perc + fl_cd.y  + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODb, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmod2b')
#Bb<-readRDS('/home/kate/CA/data/final/Bmod2b')
rm(B_OUT)

#full model remove ag_dens
B_MODc <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  gw  + ag_perc + fl_cd.y  + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODc, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmodc2')
Bc<-readRDS('/home/kate/CA/data/final/Bmodc2')
rm(B_OUT)

#full model remove ag_perc
B_MODd <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  gw + ag_dens  + fl_cd.y  + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODd, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmodd2')
Bd<-readRDS('/home/kate/CA/data/final/Bmodd2')
rm(B_OUT)

#full model remove fl_cd
B_MODe <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  gw + ag_dens + ag_perc   + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODe, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmode2')
Be<-readRDS('/home/kate/CA/data/final/Bmode2')
rm(B_OUT)

#full model remove gw
B_MODf <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  ag_dens + ag_perc + fl_cd.y  + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODf, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmodf')
Bf<-readRDS('/home/kate/CA/data/final/Bmodf')
rm(B_OUT)

#full model plus previous year landuse
B_MODg <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  gw + ag_dens + ag_perc + fl_cd.y+  factor(aglulc.y) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODg, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmodg')
Bg<-readRDS('/home/kate/CA/data/final/Bmodg')
rm(B_OUT)

#full model plus previous year landuse, no gw
B_MODh <-y_B ~ 1 + agrip_perc + agpre_perc + agapprop_perc + spi.y + agrip_perc*spi.y + agpre_perc*spi.y + agapprop_perc*spi.y + 
  ag_dens + ag_perc + fl_cd.y+  factor(aglulc.y) + factor(year) + f(HUC, model="bym", graph=H.adj, scale.model=TRUE)
B_OUT <- inla (B_MODh, family = "binomial", data=ds_B, verbose=TRUE, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE))
saveRDS(B_OUT, file='/home/kate/CA/data/final/Bmodh')
Bh<-readRDS('/home/kate/CA/data/final/Bmodh')
rm(B_OUT)



############################
#TVP Model Ouput Processing
############################
E<-readRDS('/home/kate/CA/data/E')

data<-ds
data$bins<-cut(data$spi,breaks=15)
library(ggplot2)
ggplot(data,aes(x=bins,y=px_tvp))+stat_summary(fun.y='mean',geom="point")

#extract precisions for random effects, convert to variance, and extract expected value (ie mean variance)
prec<-E$marginals.hyperpar$`Precision for the Gaussian observations`
E_eijk<-inla.emarginal(fun=function(x) 1/x, marg=prec)
prec2<-E$marginals.hyperpar$`Precision for HUC (iid component)` 
E_v00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2) 
prec2b<-E$marginals.hyperpar$`Precision for HUC (spatial component)` 
E_s00k<-inla.emarginal(fun=function(x) 1/x,marg=prec2b) 
prec3<-E$marginals.hyperpar$`Precision for POINTID` 
E_u0jk<-inla.emarginal(fun=function(x) 1/x,marg=prec3)


totMLMvarE<-E_eijk+E_v00k+E_u0jk #multilevel nonstructured variance 

#compute random intercept by watershed (level 3)

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

#plot variance of the intercept by HUC 
shp_dup<-shp
shp_dup<-merge(shp_dup, beta0jk_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
#spplot(shp_dup, "median", main = "Posterior Median of Random Intercept")#plot the  median intercept
spplot(shp_dup, c("Lower","Median","Upper"), main = "Posterior Median and 95% Credibility Interval of Random Intercept")


#compute random intercept for models with level 2 (field) and level 3 (watershed) random effects
#b0jk = b000 + v00k + u0jk

E_b000<-inla.rmarginal(1000,marg=E$marginals.fixed$`(Intercept)`) #random draw for the marginal distribution of the mean effect (intercept) 

#setup id sequence needed for random draw loop
id2<-as.data.frame(unique(ds$POINTID))#list of unique POINTIDs (level 2 units)
n2<- seq(1,length(id2$`unique(ds$POINTID)`))#count of level 2 units
id2<-cbind(n2,id2)
ds<-left_join(ds,id2, by=c("POINTID"= "unique(ds$POINTID)")) #create a column in the dataset that links the original POINTID to the simplified sequence


#build a matrix of level 2 and level 3 combined random effects
E_u0jk<-matrix(NA,1000,length(n2))#an empty matrix to hold the random effects of POINTID and HUC on the intercept term b0jk
E_U<-matrix(NA,1000,length(n2)) #empty matrix to hold the level 2 and level 3 random effects
for(i in 1:length(n2)){
  E_u0jk[ ,i]<-inla.rmarginal(1000,marg=E$marginals.random$POINTID[[i]]) #random draw from the marginal distribtuion of the level 2 RE
  n3=first(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located
  n4=last(ds$HUC[ds$n2 == i]) #level 3 unit in which the level 2 unit is located 
  E_sv00k_<-inla.rmarginal(1000,marg=(E$marginals.random$HUC[1:length(HUC)][[n3]])  #random draw from the marginal distribution of the level 3 full RE for the group in which the level 2 unit is nested
                           E_s00k<-inla.rmarginal(1000,marg= (E$marginals.random$HUC[1:length(HUC)][[n4]]) #random draw from the marginal distribution of the level 3 spatial RE for the group in which the level 2 unit is nested
                                                  E_v00k<-E_sv00k-E_s00k                       #compute the difference between the full area specific effect and the spatially structured effect to retain only iid component
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
#spplot(Eplot, c("Lower","Median","Upper"), main = "E: Posterior Median and 95% Credibility Interval of Random Intercept", cex =0.1)
plotE<- spplot(Eplot, c("Median"), key.space= list(x=0.5, y=0.95, corner=c(0,1)), main = "E: Posterior Median of Random Intercept", cex =0.2, cuts=c(-2,-1,-0.5,0.5,1,2))
plotE$legend$inside$args$key$points$cex<-c(1,1,1,1)
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
shp_dup<-merge(shp_dup, E_s00k_q, by=c("HUC"= "ID"))
colnames(shp_dup@data)<-c("HUC12", "HUC","Lower", "Median", "Upper")
spplot(shp_dup, "Median", main = "E: Posterior Median of Watershed Spatial Effect", at=c(-2,-1,-0.5,0.5,1,2))#plot the  median intercept
#spplot(shp_dup, c("Lower","Median","Upper"), main = "E: Posterior Median and 95% Credibility Interval of Watershed Spatial Effect")


#caluclate the residuals
Eres<-ds[ ,c("px_tvp","POINTID", "lat","lon", "year")]
Eres<-cbind(Eres,E$summary.fitted.values$mean,E$summary.fitted.values$sd)
colnames(Eres)<-c("px_tvp","POINTID","lat","lon","year","mean","sd")
Eres<-mutate(Eres, residuals = (px_tvp-mean)/sd)
plot(Eres$year,Eres$residuals, cex=0.5) #plot standardized residuals through time

#calculate RMSE
Ermse<-sqrt(mean((Eres$residuals*Eres$sd)^2)) 
Enrmse<-Ermse/(max(Eres$px_tvp)-min(Eres$px_tvp)) #range normalized RMSE 

#map the residuals for 1 year
coordinates(Eres)<- ~lon+lat
Eres@proj4string <- CRS('+init=epsg:3310')
Eres<-Eres[Eres@data$year==7,]
plotEres<-spplot(Eres[ ,"residuals"],cex=0.2,key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="E: Standardized \nResiduals: 2014", cuts=c(-30,-2,-1,1,2,30))
plotEres$legend$inside$args$key$points$cex<-c(1,1,1,1,1) 
print(plotEres) #plot standardized residuals through space for one year
#note


#fitted values (estimates of tvp) through time
r<-Eres
t<-r %>% group_by(year) %>% summarise(mean=mean(mean))
t2<-r %>% group_by(year) %>%summarise(sd=mean(sd))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
t$year<-t$year + 2007
plot(t$year, t$mean, ylim = c(-1,1), main="Mean Fitted TVP Values through Time", xlab="Year", ylab="Mean Fitted Values for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)
raw<-ds%>% group_by(year) %>% summarise(mean=mean(px_tvp))
raw$year<-raw$year+2007
points(raw$year, raw$mean, cex=1,col="green")


# Single Water Right-SPI Interaction Plots
moderator <-seq(-4,4) # SPI (standardized for sample)
focal<-seq(-3,3) #water right percentage (standardized)
b0<-E$summary.fixed[1,] #overall intercept
b1r<-E$summary.fixed[2,] #main effect of the focal predictor, rip in this case
b1p<-E$summary.fixed[3,] #main effect of the focal predictor, pre in this case
b1a<-E$summary.fixed[4,] #main effect of the focal predictor, approp in this case
b2<-E$summary.fixed[7,]#main effect of the moderator, spi
b3r<-E$summary.fixed[21,] #interaction effect, rip:spi
b3p<-E$summary.fixed[22,] #interaction effect, pre:spi
b3a<-E$summary.fixed[23,] #interaction effect, approp:spi

#Plot of the conditional regression of TVP on riparian rights as a fcn of SPI
plot(focal, ((b0[,1] + b2[,1]*moderator[5])+(b1r[,1] + b3r[,1]*moderator[5])*focal), type="l", ylim=c(-1,1), xlab="Percent Riparian Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Riparian Water Rights as a Function of SPI",cex.main = 1 ) #conditional regression on tvp as fcn of spi
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1r[,1] + b3r[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1r[,1] + b3r[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[3])+(b1r[,1] + b3r[,1]*moderator[3])*focal), col="red")
lines(focal, ((b0[,1] + b2[,1]*moderator[7])+(b1r[,1] + b3r[,1]*moderator[7])*focal), col="purple")
legend(0.7,1.05 , c("-2", "-1", "0", "1", "2"),title ="SPI", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("red","green", "black","blue", "purple"))

#Plot of the conditional regression of TVP on pre-1914 rights as a fcn of SPI
plot(focal, ((b0[,1] + b2[,1]*moderator[5])+(b1p[,1] + b3p[,1]*moderator[5])*focal), type="l", ylim=c(-.5,0), xlab="Percent Pre-1914 Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Pre-1914 Water Rights as a Function of SPI",cex.main = 1 ) #conditional regression on tvp as fcn of spi
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1p[,1] + b3p[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1p[,1] + b3p[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[3])+(b1p[,1] + b3p[,1]*moderator[3])*focal), col="red")
lines(focal, ((b0[,1] + b2[,1]*moderator[7])+(b1p[,1] + b3p[,1]*moderator[7])*focal), col="purple")
legend(0.7,0 , c("-2", "-1", "0", "1", "2"),title ="SPI", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("red","green", "black","blue", "purple"))

#Plot of the conditional regression of TVP on appropriative rights as a fcn of SPI
plot(focal, ((b0[,1] + b2[,1]*moderator[5])+(b1a[,1] + b3a[,1]*moderator[5])*focal), type="l", ylim=c(-1,1), xlab="Percent Appropriative Water Rights", 
     ylab="Expected Value of TVP", main="Conditional Regression of TVP on \nPercent Appropriative Water Rights as a Function of SPI",cex.main = 1 ) #conditional regression on tvp as fcn of spi
lines(focal, ((b0[,1] + b2[,1]*moderator[4])+(b1a[,1] + b3a[,1]*moderator[4])*focal), col="green")
lines(focal, ((b0[,1] + b2[,1]*moderator[6])+(b1a[,1] + b3a[,1]*moderator[6])*focal), col="blue")
lines(focal, ((b0[,1] + b2[,1]*moderator[3])+(b1a[,1] + b3a[,1]*moderator[3])*focal), col="red")
lines(focal, ((b0[,1] + b2[,1]*moderator[7])+(b1a[,1] + b3a[,1]*moderator[7])*focal), col="purple")
legend(0.7,1.05 , c("-2", "-1", "0", "1", "2"),title ="SPI", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("red","green", "black","blue", "purple"))

#plot the confidence in the simple slope of water rights on tcp as moderated by spi
#simple slope <- (b1[,1] + b3[,1]*moderator*focal[1])
mod<-as.matrix(moderator[c(2,4,5,6,8)])
b1sr<-inla.rmarginal(1000,marg=E$marginals.fixed$rip_perc)
b3sr<-inla.rmarginal(1000, marg=E$marginals.fixed$`rip_perc:spi`)
simplesloper<-matrix(NA,1000,5)
simplesloper[,1]<-(b1sr + b3sr*mod[1])
simplesloper[,2]<-(b1sr + b3sr*mod[2])
simplesloper[,3]<-(b1sr + b3sr*mod[3])
simplesloper[,4]<-(b1sr + b3sr*mod[4])
simplesloper[,5]<-(b1sr + b3sr*mod[5])
simplesloper_quartiles<-t(apply(simplesloper, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

b1sp<-inla.rmarginal(1000,marg=E$marginals.fixed$pre_perc)
b3sp<-inla.rmarginal(1000, marg=E$marginals.fixed$`pre_perc:spi`)
simpleslopep<-matrix(NA,1000,5)
simpleslopep[,1]<-(b1sp + b3sp*mod[1])
simpleslopep[,2]<-(b1sp + b3sp*mod[2])
simpleslopep[,3]<-(b1sp + b3sp*mod[3])
simpleslopep[,4]<-(b1sp + b3sp*mod[4])
simpleslopep[,5]<-(b1sp + b3sp*mod[5])
simpleslopep_quartiles<-t(apply(simpleslopep, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

b1sa<-inla.rmarginal(1000,marg=E$marginals.fixed$approp_perc)
b3sa<-inla.rmarginal(1000, marg=E$marginals.fixed$`approp_perc:spi`)
simpleslopea<-matrix(NA,1000,5)
simpleslopea[,1]<-(b1sa + b3sa*mod[1])
simpleslopea[,2]<-(b1sa + b3sa*mod[2])
simpleslopea[,3]<-(b1sa + b3sa*mod[3])
simpleslopea[,4]<-(b1sa + b3sa*mod[4])
simpleslopea[,5]<-(b1sa + b3sa*mod[5])
simpleslopea_quartiles<-t(apply(simpleslopea, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=F) #set up a right margin for the plot
plot((mod), simpleslopep_quartiles[,2], type="l", col="red", ylim=c(-0.05, 0.15), xlim=c(-3,3), ylab="Effect", xlab="SPI")
lines((mod), simpleslopep_quartiles[,1], col="red", lty=2)
lines((mod), simpleslopep_quartiles[,3], col="red", lty=2)
lines((mod), simplesloper_quartiles[,2], col="blue")
lines((mod), simplesloper_quartiles[,1], col="blue", lty=2)
lines((mod), simplesloper_quartiles[,3], col="blue", lty=2)
lines((mod), simpleslopea_quartiles[,2], col="black")
lines((mod), simpleslopea_quartiles[,1], col="black", lty=2)
lines((mod), simpleslopea_quartiles[,3], col="black", lty=2)
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=T) #allow drawing outside the plot area
legend("topright",inset=c(-0.43,0), c("Pre-1914", "Riparian","Appropriative"),title ="Type of Water Right", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("red","blue", "black"))
title( main="Effect of SPI on TVP moderated by Water Rights", cex.main=1,outer=T, line=-2) #Write the title in the top margin
par(mar=c(5, 4, 4, 2) + 0.1, xpd=F) #return plot layout settings to default

#plot the confidence in the simple slope of spi on tvp as moderated by Water Rights!!
#simple slope <- (b1[,1] + b3[,1]*moderator*focal[1])
mod<-as.matrix(focal[c(2,3,4,5,6)])
b1s<-inla.rmarginal(1000,marg=E$marginals.fixed$spi)
b3sr<-inla.rmarginal(1000, marg=E$marginals.fixed$`rip_perc:spi`)
simplesloper<-matrix(NA,1000,5)
simplesloper[,1]<-(b1s + b3sr*mod[1])
simplesloper[,2]<-(b1s + b3sr*mod[2])
simplesloper[,3]<-(b1s + b3sr*mod[3])
simplesloper[,4]<-(b1s + b3sr*mod[4])
simplesloper[,5]<-(b1s + b3sr*mod[5])
simplesloper_quartiles<-t(apply(simplesloper, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

b3sp<-inla.rmarginal(1000, marg=E$marginals.fixed$`pre_perc:spi`)
simpleslopep<-matrix(NA,1000,5)
simpleslopep[,1]<-(b1s + b3sp*mod[1])
simpleslopep[,2]<-(b1s + b3sp*mod[2])
simpleslopep[,3]<-(b1s + b3sp*mod[3])
simpleslopep[,4]<-(b1s + b3sp*mod[4])
simpleslopep[,5]<-(b1s + b3sp*mod[5])
simpleslopep_quartiles<-t(apply(simpleslopep, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))


b3sa<-inla.rmarginal(1000, marg=E$marginals.fixed$`approp_perc:spi`)
simpleslopea<-matrix(NA,1000,5)
simpleslopea[,1]<-(b1s + b3sa*mod[1])
simpleslopea[,2]<-(b1s + b3sa*mod[2])
simpleslopea[,3]<-(b1s + b3sa*mod[3])
simpleslopea[,4]<-(b1s + b3sa*mod[4])
simpleslopea[,5]<-(b1s + b3sa*mod[5])
simpleslopea_quartiles<-t(apply(simpleslopea, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=F) #set up a right margin for the plot
plot((mod), simpleslopep_quartiles[,2], type="l", col="red", ylim=c(-0.05, 0.10), xlim=c(-2,2), ylab="Effect", xlab="Percent Water Rights")
lines((mod), simpleslopep_quartiles[,1], col="red", lty=2)
lines((mod), simpleslopep_quartiles[,3], col="red", lty=2)
lines((mod), simplesloper_quartiles[,2], col="blue")
lines((mod), simplesloper_quartiles[,1], col="blue", lty=2)
lines((mod), simplesloper_quartiles[,3], col="blue", lty=2)
lines((mod), simpleslopea_quartiles[,2], col="black")
lines((mod), simpleslopea_quartiles[,1], col="black", lty=2)
lines((mod), simpleslopea_quartiles[,3], col="black", lty=2)
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=T) #allow drawing outside the plot area
legend("topright",inset=c(-0.43,0), c("Pre-1914", "Riparian","Appropriative"),title ="Type of Water Right", 
       lwd=c(2,2,2), lty=c(1,1,1), col=c("red","blue", "black"))
title( main="Effect of SPI on TVP as Moderated by Percent Water Rights", cex.main=1,outer=T, line=-2) #Write the title in the top margin
par(mar=c(5, 4, 4, 2) + 0.1, xpd=F) #return plot layout settings to default


#ALL THREE TOGETHER
b0p<-E$summary.fixed[1,] #overall intercept
b1p<-E$summary.fixed[3,] #main effect of the focal predictor, pre in this case
b2p<-E$summary.fixed[7,]#main effect of the moderator, SPI
b3p<-E$summary.fixed[22,] #interaction effect, pre:spi
b0a<-E$summary.fixed[1,] #overall intercept
b1a<-E$summary.fixed[4,] #main effect of the focal predictor, approp in this case
b2a<-E$summary.fixed[7,]#main effect of the moderator, SPI
b3a<-E$summary.fixed[23,] #interaction effect, approp:SPI
b0r<-E$summary.fixed[1,] #overall intercept
b1r<-E$summary.fixed[2,] #main effect of the focal predictor, rip in this case
b2r<-E$summary.fixed[7,]#main effect of the moderator, SPI
b3r<-E$summary.fixed[21,] #interaction effect, rip:SPI

par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=F) #set up a right margin for the plot
plot(focal, ((b0p[,1] + b2p[,1]*moderator[5] )+(b1p[,1] + b3p[,1]*moderator[5])*focal), type="l", col="red", ylim=c(-0.5,-0.1),xlim=c(-2,2), xlab="Percent Water Rights", 
     ylab="Expected Value of TVP" ) #conditional regression on tvp as fcn of time
lines(focal, ((b0p[,1] + b2p[,1]*moderator[3] )+(b1p[,1] + b3p[,1]*moderator[4])*focal), col="red", lty=2)
lines(focal, ((b0p[,1] + b2p[,1]*moderator[7])+(b1p[,1] + b3p[,1]*moderator[6])*focal), col="red", lty=3)
lines(focal, ((b0a[,1] + b2a[,1]*moderator[5])+(b1a[,1] + b3a[,1]*moderator[5])*focal), col="black")
lines(focal, ((b0a[,1] + b2a[,1]*moderator[3])+(b1a[,1] + b3a[,1]*moderator[4])*focal), col="black", lty=2)
lines(focal, ((b0a[,1] + b2a[,1]*moderator[7])+(b1a[,1] + b3a[,1]*moderator[6])*focal), col="black", lty=3)
lines(focal, ((b0r[,1] + b2r[,1]*moderator[5])+(b1r[,1] + b3r[,1]*moderator[5])*focal), col="blue")
lines(focal, ((b0r[,1] + b2r[,1]*moderator[3])+(b1r[,1] + b3r[,1]*moderator[4])*focal), col="blue", lty=2)
lines(focal, ((b0r[,1] + b2r[,1]*moderator[7])+(b1r[,1] + b3r[,1]*moderator[6])*focal), col="blue", lty=3)
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=T) #allow drawing outside the plot area
legend("topright",inset=c(-.45,0), c("Pre-1914 & Mean SPI", "Pre-1914 & -1 SD SPI","Pre-1914 & +1 SD SPI", 
                                     "Appropriative & Mean SPI","Appropriative & -1 SD SPI","Appropriative & +1 SD SPI",
                                     "Riparian & Mean SPI","Riparian & -1 SD SPI","Riparian & +1 SD SPI"),title ="Type of Water Right and \nSPI Conditions", 
       lwd=c(2,2,2,2,2,2,2,2,2), lty=c(1,2,3,1,2,3,1,2,3), col=c("red","red","red","black","black", "black","blue","blue","blue"), bty="n", cex=.75)
title( main="Conditional Regression of TVP on Percent Water Rights as a Function of SPI", cex.main=1,outer=T, line=-2) #Write the title in the top margin
par(mar=c(5, 4, 4, 2) + 0.1, xpd=F) #return plot layout settings to default


#This plot flips the x axis
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=F) #set up a right margin for the plot
plot(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[4]),xlim=c(-3,3),ylim=c(-0.45,-0.15), type="l", col="red", xlab="SPI", 
     ylab="Expected Value of TVP" ) #conditional regression on tvp as fcn of time
lines(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[3]), col="red", lty=2)
lines(moderator, ((b0p[,1] + b2p[,1]*moderator)+(b1p[,1] + b3p[,1]*moderator)*focal[5]), col="red", lty=3)
lines(moderator, ((b0a[,1] + b2a[,1]*moderator)+(b1a[,1] + b3a[,1]*moderator)*focal[4]), col="black")
lines(moderator, ((b0a[,1] + b2a[,1]*moderator)+(b1a[,1] + b3a[,1]*moderator)*focal[3]), col="black", lty=2)
lines(moderator, ((b0a[,1] + b2a[,1]*moderator)+(b1a[,1] + b3a[,1]*moderator)*focal[5]), col="black", lty=3)
lines(moderator, ((b0r[,1] + b2r[,1]*moderator)+(b1r[,1] + b3r[,1]*moderator)*focal[4]), col="blue")
lines(moderator, ((b0r[,1] + b2r[,1]*moderator)+(b1r[,1] + b3r[,1]*moderator)*focal[3]), col="blue", lty=2)
lines(moderator, ((b0r[,1] + b2r[,1]*moderator)+(b1r[,1] + b3r[,1]*moderator)*focal[5]), col="blue", lty=3)
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=T) #allow drawing outside the plot area
legend("topright", inset=c(-.45,0),  c("Mean % Pre-1914 ","-1 SD % Pre-1914 ","+1 SD % Pre-1914 ", #draw the legend in the right margin
                                       "Mean % Appropriative ", "-1 SD % Appropriative ","+1 SD % Appropriative ",
                                       "Mean % Riparian ", "-1 SD % Riparian ","+1 SD % Riparian " ),
       title ="Type and Proportion \nof Water Right", lwd=c(2,2,2,2,2,2), lty=c(1,2,3,1,2,3,1,3,3), 
       col=c("red","red","red","black","black", "black","blue","blue","blue" ), cex=0.75, bty="n")
title( main="Conditional Regression of TVP on  SPI as a Function of Percent Water Rights", cex.main=1,outer=T, line=-2) #Write the title in the top margin
par(mar=c(5, 4, 4, 2) + 0.1, xpd=F) #return plot layout settings to default


#############################################
#Barren and Fallow Model Output Processing
#############################################
#see page 140 of R-INLA book for exponentiating for interp. use antilogit for intercept and exponentiate for betas

H<-readRDS('/home/kate/CA/data/H') 

##return the antilogit of the intercept posterior marginal --> interpret as the average prob of y=1 (B&F) when all predictors at reference value or zero
Hint<-inla.tmarginal(function(x) exp(x)/(1+exp(x)),H$marginals.fixed[[1]])#build the transformed posterior
Hintquant<-inla.zmarginal(Hint) #obtain the summary stats for the intercept --> median prob of being B&F is 0.0017 or 1 out of 1,700 

##return the exponentiated posterior of the coefficients --> interpret as the incremental change in the prob of B&F given a one unit increase in predictor
#H_probincr<-inla.emarginal(exp,H$marginals.fixed$pre_perc) #this is only for the mean, recommended to use the entire posterior to avoid biases
Hrip<-inla.tmarginal(function(x) exp(x),H$marginals.fixed$rip_perc)#build the transformed posterior for riparian
Hripquant<-inla.zmarginal(Hrip)

Hpre<-inla.tmarginal(function(x) exp(x),H$marginals.fixed$pre_perc)
Hprequant<-inla.zmarginal(Hpre)

Happ<-inla.tmarginal(function(x) exp(x),H$marginals.fixed$approp_perc)
Happquant<-inla.zmarginal(Happ)

Hspi<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$spi)
Hspiquant<-inla.zmarginal(Hspi)

Hrs<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`rip_perc:spi`)
Hrsquant<-inla.zmarginal(Hrs)

Hps<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`pre_perc:spi`)
Hpsquant<-inla.zmarginal(Hps)

Has<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`approp_perc:spi`)
Hasquant<-inla.zmarginal(Has)


######################################
#Area Map
#############################

##FIGURE 1
#central valley spatail extents map
library(ggmap)
library(RColorBrewer) 

colors <- brewer.pal(9, "RdPu")

mapImage <- get_map(location = c(lon = -121, lat = 37.25),
                    color = "color",
                    source = "google",
                    maptype = "terrain",
                    zoom = 7)

cv.map <- spTransform(cvs, '+init=epsg:4326 +proj=longlat
                      +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
cv.map<-fortify(cv.map)


ggmap(mapImage, maprange=FALSE) +
  geom_polygon(aes(x = long,
                   y = lat,
                   group = group),
               data = cv.map,
               size=0.1,
               color = colors[9],
               fill = colors[6],
               alpha = 0.05) +
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
#Raw data plots and Maps
########################################

#Raw data maps

#pull out a single year from the final dataset to build coord system or will get jumbled
d3<-dsdup[dsdup$year==3,]
coordinates(d3)<-cbind(d3$lon,d3$lat)
sub<-d3


spplot(sub[ ,"AreaSqKm"], cex=0.5)
spplot(sub[ ,"wr_cnt"], cex=0.5)
spplot(sub[ ,"wr_dens"], cex=0.5)
spplot(sub[ ,"pre_cnt"], cex=0.5)
spplot(sub[ , "ag_cnt"], cex=0.5)
spplot(sub[ , "adjud_cnt"], cex=0.5)
spplot(sub[ , "spi"], cex=0.1)
spplot(sub[!is.na(sub$gw) , "gw"], cex=0.5, main='Depth to Groundwater (ft) in 2014',key.space= list(x=0.5, y=0.95, corner=c(0,1))) 
spplot(sub[ ,"diverse"])
spplot(sub[ ,"fl_cd"])
spplot(sub[ ,"adjud_perc"])
spplot(sub[!is.na(sub@data$pre_perc) ,"pre_perc"], cex=0.5)
spplot(sub[ ,"px_tvp"],cex=0.2, main='Total Vegetative Production in 2014', key.space= list(x=0.5, y=0.95, corner=c(0,1)))
spplot(sub[sub@data$px_tvp==0, ])
spplot(sub[!is.na(sub$agrip_perc) ,"agrip_perc"],cex=0.5, main='Percent Riparian in 2014', key.space= list(x=0.5, y=0.95, corner=c(0,1))) 
spplot(sub[!is.na(sub$agpre_perc) ,"agpre_perc"],cex=0.5, main='Percent Pre-1914 in 2014',key.space= list(x=0.5, y=0.95, corner=c(0,1)))
spplot(sub[!is.na(sub$agapprop_perc) ,"agapprop_perc"],cex=0.5,main='Percent Appropriative in 2014',key.space= list(x=0.5, y=0.95, corner=c(0,1)))
spplot(sub[!is.na(sub$ag_dens) ,"ag_dens"],cex=0.5)

#Landuse map with legend
z1<- spplot(sub[, "aglulc"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),
            legendEntries =c(" ","Barren & Fallow", "Grasses", "Grains", "Row Crops","Fruits and Nuts","Uncultivated Cover"),
            main="Land Use Classification in 2014")
z1$legend$inside$args$key$points$cex<-c(0,1,1,1,1,1,1)
plot(z1)

#TVP histogram
hist(ds$px_tvp, main="Histogram of All TVP Values", xlab="TVP")

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
shp<-merge(shp, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","tvp")
spplot(shp[ ,3])

#tvp mean trends through space time
st<-ds[ds$year==0, ]
st1<- st%>%group_by(HUC) %>% summarise(mean(px_tvp))
shp<-merge(shp, st1, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","tvp", "tvp0")
spplot(shp[ ,4])

st<-ds[ds$year==7, ]
st2<- st%>%group_by(HUC) %>% summarise(mean(px_tvp))
shp<-merge(shp, st2, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","tvp", "TVP_2007", "TVP_2014")
spplot(shp[ ,5])

shp@data<-mutate(shp@data, Delta_TVP=(tvp7 - tvp0))
spplot(shp[ ,6])

#Map of mean TVP by watershed for 2 years
spplot(shp,c("TVP_2007","TVP_2014"), main="Mean TVP for each Watershed")

#Maps of TVP for 2 years
sub<-ds_
sub<-sub[sub@data$year==14,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]#pull the farmland area
z1<- spplot(sub[, "px_tvp"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="TVP in 2014", cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z1$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z1)
sub<-d
sub<-sub[sub@data$year==7,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]
z2<- spplot(sub[!is.na(sub@data$px_tvp), "px_tvp"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),  main="TVP in 2007", cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z2$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z2)
print(z2, position=c(0,0,0.5,1), more=TRUE)
print(z1, position=c(0.5,0,1,1))

#Groundwater maps for 2 years
sub<-d
sub<-sub[sub@data$year==14,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]#pull the farmland area
z1<- spplot(sub[, "gw"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)), main="TVP in 2014")#, cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z1$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z1)
sub<-d
sub<-sub[sub@data$year==7,] #pull out one year at a time
sub<-sub[sub@data$flflag==1,]
z2<- spplot(sub[, "gw"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),  main="TVP in 2007")#, cuts=c(0, 0.169, 0.333, 0.497, 0.662, 0.826,1.337))
z2$legend$inside$args$key$points$cex<-c(1,1,1,1,1)
#plot(z2)
print(z2, position=c(0,0,0.5,1), more=TRUE)
print(z1, position=c(0.5,0,1,1))


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
shp2<-merge(shp2, s, by="HUC")
colnames(shp2@data)<-c("HUC12", "HUC","spi")
spplot(shp2[ ,3])

#spi mean trends through space time
st<-ds[ds$year==0, ]
st1<- st%>%group_by(HUC) %>% summarise(mean(spi))
shp2<-merge(shp2, st1, by="HUC")
colnames(shp2@data)<-c("HUC12", "HUC","spi", "spi2007")
spplot(shp2[ ,4])

st<-ds[ds$year==7, ]
st2<- st%>%group_by(HUC) %>% summarise(mean(spi))
shp2<-merge(shp2, st2, by="HUC")
colnames(shp2@data)<-c("HUC12", "HUC","SPI", "SPI_2007", "SPI_2014")

#Maps of SPI by watershed for 2 years
spplot(shp2, c("SPI_2007","SPI_2014"), main ="Mean SPI for each Watershed")

#Ag water rights density trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(ag_dens))
s2<-ds %>% group_by(HUC) %>%summarise(sd(ag_dens))
shp<-merge(shp, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","Water_Rights_Density")
spplot(shp[ ,"Water_Rights_Density"])

#Ag water rights density mean trend through time
t<-ds %>% group_by(year) %>% summarise(mean=mean(ag_dens, na.rm=T))
t2<-ds %>% group_by(year) %>%summarise(sd=sd(ag_dens, na.rm=T))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim = c(-7,4), main="Mean WR Density through Time", xlab="Years since 2007", ylab="Mean SPI for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

plot(ds$year, ds$ag_dens)

#Riparian trend through space
s<-ds %>% group_by(HUC) %>%summarise(mean(agrip_perc, na.rm=T))
s2<-ds %>% group_by(HUC) %>%summarise(sd(agrip_perc, na.rm=T))
shp<-merge(shp, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC","Percent_Riparian")
spplot(shp[ ,"Percent_Riparian"])

t<-d3@data %>% group_by(year) %>% summarise(mean=mean(agrip_perc, na.rm=T))
t2<-d3@data %>% group_by(year) %>%summarise(sd=sd(agrip_perc, na.rm=T))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim=c(0,100), main="Mean Rip Perc through Time", xlab="Years since 2007", ylab="Mean SPI for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

plot(ds$year, ds$agrip_perc)

#Appropriative trend through space
shp<-shpdup
s<-ds %>% group_by(HUC) %>%summarise(mean(agapprop_perc, na.rm=T))
s2<-ds %>% group_by(HUC) %>%summarise(sd(agapprop_perc, na.rm=T))
shp<-merge(shp, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC", "Percent_Post_1914_Appropriative")
#spplot(shp[ ,3], main="Percent Post-1914 Appropriative")

t<-d3@data %>% group_by(year) %>% summarise(mean=mean(agapprop_perc, na.rm=T))
t2<-d3@data %>% group_by(year) %>%summarise(sd=sd(agapprop_perc, na.rm=T))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim=c(0,100), main="Mean Approp Perc through Time", xlab="Years since 2007", ylab="Mean SPI for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

plot(ds$year, ds$agapprop_perc)

#Pre-1914 trend through space
shp<-shpdup
s<-ds %>% group_by(HUC) %>%summarise(mean(agpre_perc))
s2<-ds %>% group_by(HUC) %>%summarise(sd(agpre_perc))
shp<-merge(shp, s, by="HUC")
colnames(shp@data)<-c("HUC12", "HUC", "Percent_Pre_1914")
spplot(shp[ ,3], main="Percent Pre-1914 Appropriative")

t<-d3@data %>% group_by(year) %>% summarise(mean=mean(agpre_perc, na.rm=T))
t2<-d3@data %>% group_by(year) %>%summarise(sd=sd(agpre_perc, na.rm=T))
t<-mutate(t,upper=t$mean+t2$sd, lower=t$mean-t2$sd)
plot(t$year, t$mean, ylim=c(0,100), main="Mean Pre Perc through Time", xlab="Years since 2007", ylab="Mean SPI for the Central Valley", type="l")
lines(t$year,t$upper, col="red", lty=2)
lines(t$year,t$lower, col="red", lty=2)

plot(ds$year, ds$agpre_perc)


##Map of land use
sub<-d
sub<-sub[sub@data$flflag==1,]
sub<-sub[sub@data$year==14,] #pull out one year at a time
z1<- spplot(sub[, "aglulc"], cex=0.1, key.space= list(x=0.5, y=0.95, corner=c(0,1)),
            legendEntries =c(" ","Barren & Fallow", "Grasses", "Grains", "Row Crops","Fruits and Nuts","Uncultivated Cover"),
            main="Land Use Classification in 2014")
z1$legend$inside$args$key$points$cex<-c(0,1,1,1,1,1,1)
plot(z1)

#Basic stats
mean(d$ag_dens, na.rm=T)
sd(d$px_tvp, na.rm=T)
min(ds$gw, na.rm=T)
max(dsdup$ag_dens, na.rm=T)
