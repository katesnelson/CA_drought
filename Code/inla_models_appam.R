library(memisc)
library(dplyr)
library (gstat)
library(sp)
library(nlme)
library(rgdal)
library(INLA)
library(spdep)
library(raster)

setwd('/home/kate/R/x86_64-pc-linux-gnu-library/3.3/INLA/')
##############################################
#Read in datasets to use for inla models
##############################################
#full cleaned pixel dataset of central valley as spatial points data frame
#d <- readRDS('/data/emily/WF/kate/final_data/cv.dat.clean.rds')

###  full dataset with lat lon columns and updated wr info --> check for cleanliness 
d<- readRDS('/home/kate/CA/data/cv.dat2.rds') 
#plot(d)


#full central valley watershed as spatial polygons data frame
cv <- shapefile('/data/emily/WF/kate/shp/cv_huc.shp')
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
m=mean(ds$px_tvp)
stdev=sd(ds$px_tvp)
y1=(ds$px_tvp - m)/stdev #standardized dependent variable

ds$year <- ds$year - 7 #center time at 2007 = 0
ds <- mutate(ds, year1=year) #add duplicate year grouping index
ds<-mutate(ds,POINTID1=POINTID)#add duplicate pixel grouping index
ds <- mutate(ds, HUC1=HUC,HUC2=HUC, HUC3=HUC, HUC4=HUC, HUC5=HUC, HUC6=HUC, HUC7=HUC, HUC8=HUC,HUC9=HUC, 
             HUC10=HUC, HUC11=HUC, HUC12=HUC, HUC13=HUC, HUC14=HUC, HUC15=HUC, HUC16=HUC, HUC17=HUC) #add duplicates of HUC index for random slope grouping

s.index1 <- inla.spde.make.index(name="spatial.field1", n.spde=spde$n.spde)

##################################################
#Models with no Predictors and basic error structure
##################################################
#adding control.predictor(compute=T) caluclates the marginal posterior distribution for fitted data values of the linear predictor (good for maps of fit at points)but takes a long time to run

##Null 2 level MLM -huc grouped ##  takes approx 5 min to 15min to run --> essentially aggregates pixel variation
null <-y ~ 1 + f(HUC, model="iid")
output.null <- inla (null, family = "gaussian", data=ds) #takes very long time to run when add control.compute = list(dic = T) or waic=T
saveRDS(output.null, file='/home/kate/CA/data/output.null.huc')
nh<-readRDS('/home/kate/CA/data/output.null.huc2')
# null.summary <- summary (output.null) #
# ICC.null <- null.summary$hyperpar$mean[2]/(null.summary$hyperpar$mean[1]+null.summary$hyperpar$mean[2])#0.48 - lots of between huc error

##Null space-time model ##  
null3 <-y ~ 1 + f(HUC, year, model="iid") # space-time error structure
output.null3 <- inla (null3, family = "gaussian", data=ds) 
                     # control.predictor =list( compute=TRUE))
saveRDS(output.null3, file='/home/kate/CA/data/output.null.st')
nst<-readRDS('/home/kate/CA/data/output.null.st')
# null.summary3 <- summary (output.null3) #
# ICC.null3 <- null.summary3$hyperpar$mean[2]/(null.summary3$hyperpar$mean[1]+null.summary3$hyperpar$mean[2]) #0.95 - most error is between years

##Null 3 level MLM ## 3 level model for instances (time) nested hucs and pixels nested in hucs
null_3lvl  <-y ~ 1 + f(HUC, model="iid") + f(HUC1,POINTID, model = "iid")+ f(HUC2,year, model ="iid")
output.null.3lvl <- inla (null_3lvl, family = "gaussian", data=ds, verbose=TRUE)
saveRDS(output.null.3lvl, file='/home/kate/CA/data/output.null.3lvl')
n3<-readRDS('/home/kate/CA/data/output.null.3lvl')

##Null 3 level MLM ## 3 level model for instances (time) nested in pixels nested in hucs, use to get 3 level ICC --1.5 hr run time
null_3lvl2  <-y ~ 1 + f(HUC, model="iid") + f(POINTID, model = "iid")
output.null.3lvl2 <- inla (null_3lvl2, family = "gaussian", data=ds, verbose=TRUE)
saveRDS(output.null.3lvl2, file='/home/kate/CA/data/output.null.3lvl3')
n33<-readRDS('/home/kate/CA/data/output.null.3lvl3')

##################################################################
#Models with no predictors and complex error structures
#############################################################

##HUC spatial model, level 3## for 1 year takes ~5 min, for all years takes ~1hr
huc <-y ~ f(HUC, model="besag", graph=H.adj)
output.huc <- inla (huc, family ="gaussian", data=ds, verbose=TRUE)
saveRDS(output.huc, file='/home/kate/CA/fidata/output.huc.besag')
besag<-readRDS('/home/kate/CA/data/output.huc.besag')

##HUC spatial model 2, level 3##  for all years takes ~8min
huc <-y ~ f(HUC, model="bym", graph=H.adj) 
output.huc <- inla (huc, family ="gaussian", data=ds)
saveRDS(output.huc, file='/home/kate/CA/data/output.huc.bym')
bym<-readRDS('/home/kate/CA/data/output.huc.bym')
# huc.summary2 <-summary(output.huc2) #increased marginal log-likelihood (355734.07 ) over huc null model (355047.98 )and besag model

##HUC spatial model with scaling, level 3##  for all years takes ~8min
huc2 <-y ~ f(HUC, model="bym", graph=H.adj,scale.model=TRUE) 
output.huc2 <- inla (huc2, family ="gaussian", data=ds)
saveRDS(output.huc2, file='/home/kate/CA/data/output.huc.bym2')
bym2<-readRDS('/home/kate/CA/data/output.huc.bym2') #provides slightly higher log-likelihood over nonscaled model


## Spatio-temporal models - parametric for time, levels 1 & 3 ##
temphuc3<-y~ 1 + f(HUC, model="bym", graph = H.adj,scale.model=TRUE) + f(HUC1, year, model="iid", constr=TRUE)
output.temphuc3 <- inla (temphuc3,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.temphuc3, file='/home/kate/CA/data/output.temphuc3')
temphuc3<-readRDS('/home/kate/CA/data/output.temphuc3')

temphuc6<-y~ 1 + f(HUC, model="bym", graph = H.adj,scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE)
output.temphuc6 <- inla (temphuc6,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.temphuc6, file='/home/kate/CA/data/output.temphuc6')
temphuc6<-readRDS('/home/kate/CA/data/output.temphuc6')

#add a time fixed effect
temphuc4<-y~ 1 + f(HUC, model="bym", graph = H.adj,scale.model=TRUE) + f(HUC1, year, model="iid", constr=TRUE) +year
output.temphuc4 <- inla (temphuc4,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.temphuc4, file='/home/kate/CA/data/output.temphuc4')
temphuc4<-readRDS('/home/kate/CA/data/output.temphuc4')

temphuc5<-y~ 1 + f(HUC, model="bym", graph = H.adj,scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) +year
output.temphuc5 <- inla (temphuc5,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.temphuc5, file='/home/kate/CA/data/output.temphuc5')
temphuc5<-readRDS('/home/kate/CA/data/output.temphuc5')

## Temporal and Pixel spatio-temporal model, levels 1 & 2 ##
temppix <-y ~ -1 + intercept + f(spatial.field, model = spde) + f(spatial.field1, year, model="ar1") 
temppix.stack.est <- inla.stack(data =list(y=y), 
                                A=list(A.est,1),
                                effects=list(c(s.index,s.index1, list(intercept=1)), 
                                             list(year=ds$year)),
                                tag="est" )
output.temppix <- inla (temppix,  data=inla.stack.data(temppix.stack.est, spde=spde), 
                        control.predictor =list(A=inla.stack.A(temppix.stack.est)), verbose=TRUE)
saveRDS(output.temppix, file='/home/kate/CA/data/output.temppix') #fit is much better for the huc space-time interaction model


## Temporal and Pixel and HUC spatio-temporal model, levels 1, 2 & 3 ##
#this model accounts for spatial effects a the pixel level usning an spde model, spatial effects at the huc level using a CAR + iid model, inlcudes a time fixed effect, and allows the spatial effect of hucs to vary with time according to ar1
temppixhuc<-y ~ -1 + intercept + f(spatial.field, model = spde) + f(HUC, model="bym", graph = H.adj, scale.model=TRUE)  + f(HUC1, year, model="ar1", constr=TRUE) + year
temppixhuc.stack.est <- inla.stack(data =list(y=y), 
                                   A=list(A.est,1,1,1),
                                   effects=list(c(s.index, list(intercept=1)), 
                                                list(year=ds$year),list(HUC=ds$HUC), list(HUC1=ds$HUC1)),
                                   tag="est")
output.temppixhuc <- inla (temppixhuc,  data=inla.stack.data(temppixhuc.stack.est, spde=spde), 
                           control.predictor =list(A=inla.stack.A(temppixhuc.stack.est)))
saveRDS(output.temppixhuc, file='/home/kate/CA/data/output.temppixhuc') #1.5 hr run time, better fit than without level 2 effects

#this model accounts for pixel level variance, spatial effects at the huc level using a CAR + iid model a time fixed effect and allows the spatial effect of hucs to vary with time according to ar1
temppixhuc2<-y ~ -1 + intercept + f(POINTID, model = "iid") + f(HUC, model="bym", graph = H.adj, scale.model=TRUE)  + f(HUC1, year, model="ar1") + year
output.temppixhuc2 <- inla (temppixhuc2,  data=ds, family="gaussian", verbose=TRUE)
saveRDS(output.temppixhuc2, file='/home/kate/CA/data/output.temppixhuc2')


############################################################
#Models with time and space error structures and predictors -no interaction and no random slopes
############################################################

# model 1- controls and time only
mlm1<-y~ 1 + year  + gw + spi + factor(aglulc) + fl_cd+ 
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) 
output.mlm1 <- inla (mlm1,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.mlm1, file='/home/kate/CA/data/output.mlm1')
mlm1<-readRDS('/home/kate/CA/data/output.mlm1')

#model 2- controls and time and wr_dens
mlm2<-y~ 1 + year + wr_dens + gw + spi + factor(aglulc) + fl_cd+ 
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) 
output.mlm2 <- inla (mlm2,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.mlm2, file='/home/kate/CA/data/output.mlm2')
mlm2<-readRDS('/home/kate/CA/data/output.mlm2')

# model 5- controls and time, wr_dens, rip_perc, approp_perc, pre_perc
mlm5<-y~ 1 + year + wr_dens + rip_perc + approp_perc +pre_perc + gw + spi + factor(aglulc) + fl_cd+ 
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) 
output.mlm5 <- inla (mlm5,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.mlm5, file='/home/kate/CA/data/output.mlm5')
mlm5<-readRDS('/home/kate/CA/data/output.mlm5')


##############################################################################
#MLM with space-time error and predictors with  interactions with time - no random slopes
############################################################################
intmlm<-y~ 1 +  wr_dens + rip_perc + approp_perc +pre_perc + year*wr_dens +year*rip_perc +year*approp_perc + year*pre_perc + #predictors and interactions
  year + gw + spi + factor(aglulc) + fl_cd + #controls- should we have our controls interact with time?
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) #error structure
output.intmlm <- inla (intmlm,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.intmlm, file='/home/kate/CA/data/output.intmlm')
intmlm<-readRDS('/home/kate/CA/data/output.intmlm')

intmlm2<-y~ 1 +  wr_dens + rip_perc + approp_perc +pre_perc + year*wr_dens +year*rip_perc +year*approp_perc + year*pre_perc + #predictors and interactions
  year + gw + spi + factor(aglulc) + fl_cd + gw*year + fl_cd*year + spi*year + #factor(aglulc*year)+ #controls and control interactions with time 
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) #error structure
output.intmlm2 <- inla (intmlm2,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.intmlm2, file='/home/kate/CA/data/output.intmlm2')
intmlm2<-readRDS('/home/kate/CA/data/output.intmlm2')


##############################################################################
#MLM with spcae-tiem error, predictors and interactions, and random slopes
#############################################################################
rsmlm<-y~ 1 +  wr_dens + rip_perc + approp_perc +pre_perc + year*wr_dens +year*rip_perc +year*approp_perc + year*pre_perc+ #predictors and interactions
  year + gw + spi + factor(aglulc) + fl_cd+ #controls- 
  f(HUC2, wr_dens, model="iid") + f(HUC3, rip_perc, model="iid") + f(HUC4, approp_perc, model="iid") + f(HUC5, pre_perc, model="iid"), #allow predictors to have random slope by HUC
f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) #error structure
output.rsmlm <- inla (rsmlm,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.rsmlm, file='/home/kate/CA/data/output.rsmlm')
rsmlm<-readRDS('/home/kate/CA/data/output.rsmlm')

rsmlm2<-y~ 1 +  wr_dens + rip_perc + approp_perc +pre_perc + year*wr_dens +year*rip_perc +year*approp_perc + year*pre_perc+ #predictors and interactions
  year + gw + spi + factor(aglulc) + fl_cd + gw*year + fl_cd*year + spi*year + #factor(aglulc*year), #controls and control interactions with time 
  f(HUC2, wr_dens, model="iid") + f(HUC3, rip_perc, model="iid") + f(HUC4, approp_perc, model="iid") + f(HUC5, pre_perc, model="iid")+ #allow predictors to have random slope (slope varies by HUC)
  f(HUC6, gw, model="iid") + f(HUC7, fl_cd, model="iid") + f(HUC8, spi, model="iid")+ f(HUC9, factor(aglulc), model="iid")+ #allow controls to have random slope (slope varies by HUC)
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) #error structure
output.rsmlm2 <- inla (rsmlm2,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.rsmlm2, file='/home/kate/CA/data/output.rsmlm2')
rsmlm2<-readRDS('/home/kate/CA/data/output.rsmlm2')

rsmlm3<-y~ 1 +  wr_dens + rip_perc + approp_perc +pre_perc + year*wr_dens +year*rip_perc +year*approp_perc + year*pre_perc+ #predictors and interactions
  year + gw + spi + factor(aglulc) + fl_cd + gw*year + fl_cd*year + spi*year + #factor(aglulc*year), #controls and control interactions with time 
  f(HUC2, wr_dens, model="iid") + f(HUC3, rip_perc, model="iid") + f(HUC4, approp_perc, model="iid") + f(HUC5, pre_perc, model="iid")+ #allow predictors to have random slope (slope varies by HUC)
  f(HUC6, gw, model="iid") + f(HUC7, fl_cd, model="iid") + f(HUC8, spi, model="iid")+ f(HUC9, factor(aglulc), model="iid")+#allow controls to have random slope (slope varies by HUC)
  f(HUC10, year*wr_dens, model="iid") + f(HUC11, year*rip_perc, model="iid") + f(HUC12, year*approp_perc, model="iid") + f(HUC13, year*pre_perc, model="iid")+ #allow the predictor interactions with time to vary by HUC
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) #error structure
output.rsmlm3 <- inla (rsmlm3,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.rsmlm3, file='/home/kate/CA/data/output.rsmlm3')
rsmlm3<-readRDS('/home/kate/CA/data/output.rsmlm3')

rsmlm4<-y~ 1 +  wr_dens + rip_perc + approp_perc +pre_perc + year*wr_dens +year*rip_perc +year*approp_perc + year*pre_perc, #predictors and interactions
year + gw + spi + factor(aglulc) + fl_cd + gw*year + fl_cd*year + spi*year + factor(aglulc)*year, #controls and control interactions with time 
f(HUC2, wr_dens, model="iid") + f(HUC3, rip_perc, model="iid") + f(HUC4, approp_perc, model="iid") + f(HUC5, pre_perc, model="iid")+ #allow predictors to have random slope (slope varies by HUC)
  f(HUC6, gw, model="iid") + f(HUC7, fl_cd, model="iid") + f(HUC8, spi, model="iid")+ f(HUC9, factor(aglulc), model="iid")+#allow controls to have random slope (slope varies by HUC)
  f(HUC10, year*wr_dens, model="iid") + f(HUC11, year*rip_perc, model="iid") + f(HUC12, year*approp_perc, model="iid") + f(HUC13, year*pre_perc, model="iid")+ #allow the predictor interactions with time to vary by HUC
  f(HUC14, gw*year, model="iid") + f(HUC15, fl_cd*year, model="iid") + f(HUC16, spi*year, model="iid") + #f(HUC17,factor(aglulc)*year, model="iid")+ # allow control interactions with time to vary by HUC
  f(HUC, model="bym", graph = H.adj, scale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) #error structure
output.rsmlm4 <- inla (rsmlm4,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.rsmlm4, file='/home/kate/CA/data/output.rsmlm4')
rsmlm4<-readRDS('/home/kate/CA/data/output.rsmlm4')

########################################
#Test method for letting predictors have non linear trends
#################################################

test<-y~ 1 + year + wr_dens + f(HUC2, rip_perc, model="rw2") + approp_perc +pre_perc + gw + spi + factor(aglulc) + fl_cd + 
  f(HUC, model="bym", graph = H.adj, sscale.model=TRUE) + f(HUC1, year, model="ar1", constr=TRUE) 
output.test <- inla (test,  family = "gaussian", data=ds, verbose=TRUE)#, 
#control.predictor =list( compute=TRUE))
saveRDS(output.test, file='/home/kate/CA/data/output.test')
test<-readRDS('/home/kate/CA/data/output.test')

###############################
#Plots
#############################


#linear trend of time
x <- seq(1,8) # Years
plot(x,intmlm5$summary.fixed[2,1]*x, type="l", main="",xlab="t",ylab=expression(beta*t))
lines(intmlm5$summary.fixed[2,3]*x,lty=2)
lines(intmlm5$summary.fixed[2,5]*x,lty=2)

#LINERA INTERACTION OF WR_DENS AND TIME --> INT NOT SIG
wmax<-max(ds$wr_dens)
wmin<-min(ds$wr_dens)+0.001
wmean<-mean(ds$wr_dens)
plot(x,intmlm5$summary.fixed[15,1]*x*wmax , type="l", main="",xlab="year",ylab=expression(beta*t),ylim=c(0.00001,0.003))
lines(intmlm5$summary.fixed[15,1]*x*wmin,lty=2, col="red")
lines(intmlm5$summary.fixed[15,1]*x*wmean,lty=2, col="blue")

#linear interaction of rip_perc and time
plot(x,intmlm5$summary.fixed[16,1]*x*1 + intmlm5$summary.fixed[3,1]*1 , type="l", main="",xlab="year",ylab=expression(effect of % rip on tvp))#,ylim=c(0.00001,0.003))
lines(intmlm5$summary.fixed[16,1]*x*.1 + intmlm5$summary.fixed[3,1]*0.1,lty=2, col="red")
lines(intmlm5$summary.fixed[16,1]*x*.5 + intmlm5$summary.fixed[3,1]*0.5,lty=2, col="blue")


