library(memisc)
library(dplyr)
require(ggplot2)
require(rjags)
require(ggmcmc)
require(BEST)
require(foreign)
require(shinystan)

dsf <- tbl_df(read.csv("C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Kate_Emily\\Data\\pt_data_kn.txt", stringsAsFactors = FALSE))
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
df <- ds[ds$farm_flag == 1,]
reomve(ds)
remove(dsf)

#additional variables
df$Rip <- df$COUNT_Rip/df$fmmp_area
df$P1914 <- df$COUNT_Pre1914/df$fmmp_area


#data setup
y <- df$delta_tvp
n <- length(y)

dc.names <- as.vector(df$HUC12)
uq <- unique(dc.names)
n.dc <- length(uq)
dc <- rep(NA, n.dc)
for (i in 1:n.dc){
  dc[dc.names == uq[i]] <- i
  sample.size <- as.vector(table(dc))
}

model_string <- "model{
for (i in 1:n){
  y[i] ~ dnorm(mu[i], sigma.y)
  mu[i] <- a[dc[i]] + b1*x1[dc[i]] + b2*x2[i] + b3*x3[i] + b4*x4[i] + b5*x5[i] + b6*x1[i]*x2[i]
}

b1 ~ dt(0,.1, 1) 
b2 ~ dt(0,.1, 1) 
b3 ~ dt(0,.1, 1) 
b4 ~ dt(0,.1, 1) 
b5 ~ dt(0,.1, 1) 
b6 ~ dt(0,.1, 1) 

for (j in 1:n.dc){
a[j] ~ dnorm(mu.a, tau.a)  
}

mu.a ~ dt(0, .1, 1)
sigma.y ~ dgamma(0.01, 0.01)
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

