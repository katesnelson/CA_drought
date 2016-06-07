library(R2WinBUGS)
library(dplyr)
library(maptools)
library(sp)
library(rgdal)


dsf <- tbl_df(read.csv("C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Kate_Emily\\Data\\pt_data_kn.txt", stringsAsFactors = FALSE))
#ds = dsf[dsf$farm_flag==1,]
dsf = transform(dsf, 
               tvp_09 = as.numeric(tvp_09),
               tvp_14 = as.numeric(tvp_14),
               delta_tvp = as.numeric(delta_tvp),
               AreaAcres = as.numeric(AreaAcres))
ds[ds==999.99] <- NA
dsf$Rip <- df$COUNT_Rip/df$fmmp_area
dsf$P1914 <- df$COUNT_Pre1914/df$fmmp_area
y <- dsf$delta_tvp
n <- length(y)
dc.names <- as.vector(dsf$HUC12)
uq <- unique(dc.names)
n.dc <- length(uq)
dc <- rep(NA, n.dc)
for (i in 1:n.dc){
  dc[dc.names == uq[i]] <- i
  sample.size <- as.vector(table(dc))
}

#sampling grid to create adjacency matrix
row_num = 1155
col_num = 1350
minx <- min(dsf$lon)
maxx <- max(dsf$lon)
miny <- min(dsf$lat)
maxy <- max(dsf$lat)
lon <- seq(from = minx, to = maxx, by = (abs(maxx - minx)/row_num))
lat <- seq(from = miny, to = maxy, by = (abs(maxy - miny)/col_num))
xy <- expand.grid(x = lon, y = lat)
grid.pts <- SpatialPointsDataFrame(coords = xy, data = xy, proj4string = CRS("+init=epsg:3310") )
gridded(grid.pts) <- TRUE
gridsp <- as (grid.pts, "SpatialPolygons") #this is taking forever
grid <- SpatialPolygonsDataFrame(grid.pts, data = data.frame(id=row.names(gridsp), 

#I need to import this file as a SpatialPolygonsDF, then I can build matrix
fn <- "C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Kate_Emily\\Data\\pt_final.shp"                                                    
readShapePoly(fn, IDvar = NULL, proj4string = CRS("+init=epsg:3310"),
              verbose = F, repair = F, force_ring = F, delete_null_obj = F,
              retrieve_ABS_null = F)


#spedp to compute adjacency matrix
library(spdep)
shape_nb <- poly2nb(grid.pts)
NumCells <- length(shape_nb)
num = sapply(shape_nb, length)
adj = unlist(shape_nb)
sumNumNeigh = length(unlist(shape_nb))

model_string <- "model{

for (i in 1:n){

  y[i] ~ dnorm(mu[i], sigma.y)
  mu[i] <- rr[i] * e[i]
  rr[i] <- a[dc[i]] + b5*x5[i] + v[i] + h[i]

          #h is car
          # v is random effects
          #r is residual

  r[i] <- (y[i] - mu[i])
  v[i] ~ dnorm(0, tau.v) #prior on random effects
}

#l1 priors
b5 ~ dt(0, .001, 1) 
h[1:i] ~car.normal(adj[], weights[], num[], tau.h) #prior on car
tau.v ~ dgamma(0.001, 0.001)


for (j in 1:n.dc){
a[j] ~ dnorm(mu.a, tau.a)  
}
#l2 priors
mu.a ~ dt(0, 1, 1)
sigma.y ~ dgamma(0.001, 0.001)
tau.a <- pow(sigma.a , -2)  
sigma.a ~ dunif(0, 100)  
}"

#initialize variables
inits <- function(chain) {
  list (a=rnorm(n.dc), b1 = rnorm(1), b2 = rnorm(1),
        b3 = rnorm(1), b4 = rnorm(1), b5 = rnorm(1),
        b6 = rnorm(1), sigma.y = runif(1), mu.a = runif(1), sigma.a = runif(1))}

#create dataframe
#B1 = delta_lc; B2 = diverse; B3 = WR_density; B4 = P1914; B5 = GW_dnsty15; B6 = delta_lc*diverse 
data <- list(n = n, n.dc = n.dc, y = y, dc = dc, 
             x1 = df$delta_lc, x2 = df$diverse, x3 = df$WR_density, x4 = df$P1914, x5 = df$GW_dnsty15)

#tell JAGS parameters to report back
parameters <- c("a", "b1", "b2", "b3", "b4", "b5", "b6", "mu.a", "sigma.a", "sigma.y")

#compile jags model
model <- jags.model(textConnection(model_string),
                    data = data, 
                    inits = inits,
                    n.chains = 3,
                    n.adapt = 200)

#take 2000 random samples of each of the 3 chains
update(model, n.iter = 20000, n.thin = 10, n.burnin = 5000)
model_outcome <- coda.samples(model, variable.names = parameters, n.iter = 20000)
my_sso <- as.shinystan(model_outcome)
my_sso <- launch_shinystan(my_sso)

#non cauchy priors had triple bumps in prior estimations

