
library(INLA)
library(ggplot2)
library(dplyr)

setwd('C:')

#This script build plots from INLA model output as seen in Nelson, K. S., & Burchfield, E. K. Effects of the Structure of Water Rights on Agricultural Production 
#During Drought: A Spatiotemporal Analysis of California's Central Valley. Water Resources Research.DOI: 10.1002/2017WR020666

E<-readRDS('C:/tvpmodelname')
summary(E)

H<-readRDS('C:/logisticmodelname')
summary(H)

#caterpillar plot
esf<-as.data.frame(E$summary.fixed)
colnames(esf)<-c("mean","sd","Lower","Median","Upper","mode","kld")
esf<-esf[c(2:5,21:23),]
names<-row.names(esf)
names<-c("Percent Riparian","Percent Pre-1914","Percent Appropriative","SPI","Percent Riparian*SPI","Percent Pre-1914*SPI","Percent Appropriative*SPI")
g<- ggplot(esf, aes(x= names, y=Median)) #,greek = F, thick_ci = c(0.17, 0.83))) +
g + geom_pointrange(aes(ymin=Lower, ymax=Upper),color="grey", size = 1, fatten = 0.1 ) + 
  geom_pointrange(aes(ymin=Median - sd, ymax=Median + sd),color="dark gray", size = 2, fatten = 0.1 )+
  geom_point(fill="black") + 
  scale_x_discrete(limits=c("Percent Appropriative*SPI","Percent Pre-1914*SPI","Percent Riparian*SPI","SPI","Percent Appropriative","Percent Pre-1914","Percent Riparian")) +
 coord_flip() + labs(y = "Estimated Effect on TVP", x = "Predictor") + theme_bw(base_size = 11)


#plot the confidence in the simple slope of spi on tvp as moderated by Water Rights!!
#simple slope <- (b1[,1] + b3[,1]*moderator*focal[1])
moderator <-seq(-4,4) # SPI (standardized for sample)
focal<-seq(-3,3) #water right percentage (standardized)
b0<-E$summary.fixed[1,] #overall intercept
b1r<-E$summary.fixed[2,] #main effect of  rip 
b1p<-E$summary.fixed[3,] #main effect of  pre 
b1a<-E$summary.fixed[4,] #main effect of  approp 
b2<-E$summary.fixed[5,]#main effect  spi
b3r<-E$summary.fixed[21,] #interaction effect, rip:spi
b3p<-E$summary.fixed[22,] #interaction effect, pre:spi
b3a<-E$summary.fixed[23,] #interaction effect, approp:spi
mod<-as.matrix(focal[c(2,3,4,5,6)]) #flip moderator and focal as should be thinking of water rights as moderating effect of drought

b1s<-inla.rmarginal(1000,marg=E$marginals.fixed$spi)
b3sr<-inla.rmarginal(1000, marg=E$marginals.fixed$`agrip_perc:spi`)
simplesloper<-matrix(NA,1000,5)
simplesloper[,1]<-(b1s + b3sr*mod[1])
simplesloper[,2]<-(b1s + b3sr*mod[2])
simplesloper[,3]<-(b1s + b3sr*mod[3])
simplesloper[,4]<-(b1s + b3sr*mod[4])
simplesloper[,5]<-(b1s + b3sr*mod[5])
simplesloper_quartiles<-t(apply(simplesloper, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

b3sp<-inla.rmarginal(1000, marg=E$marginals.fixed$`agpre_perc:spi`)
simpleslopep<-matrix(NA,1000,5)
simpleslopep[,1]<-(b1s + b3sp*mod[1])
simpleslopep[,2]<-(b1s + b3sp*mod[2])
simpleslopep[,3]<-(b1s + b3sp*mod[3])
simpleslopep[,4]<-(b1s + b3sp*mod[4])
simpleslopep[,5]<-(b1s + b3sp*mod[5])
simpleslopep_quartiles<-t(apply(simpleslopep, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

b3sa<-inla.rmarginal(1000, marg=E$marginals.fixed$`agapprop_perc:spi`)
simpleslopea<-matrix(NA,1000,5)
simpleslopea[,1]<-(b1s + b3sa*mod[1])
simpleslopea[,2]<-(b1s + b3sa*mod[2])
simpleslopea[,3]<-(b1s + b3sa*mod[3])
simpleslopea[,4]<-(b1s + b3sa*mod[4])
simpleslopea[,5]<-(b1s + b3sa*mod[5])
simpleslopea_quartiles<-t(apply(simpleslopea, MARGIN=2, function(x) quantile(x, probs=c(0.025,0.5,0.975))))

simpleslopea_q<-as.data.frame(simpleslopea_quartiles)
simplesloper_q<-as.data.frame(simplesloper_quartiles)
simpleslopep_q<-as.data.frame(simpleslopep_quartiles)

#and finally create the plot!
ggplot(simpleslopea_q, aes(x=mod)) +
  geom_line(aes(y=simpleslopea_q$`50%`),colour="black", lwd=1) +
  geom_ribbon(aes(ymin=simpleslopea_q$`2.5%`, ymax=simpleslopea_q$`97.5%`),alpha=0.2) +
  geom_line(aes(y=simplesloper_q$`50%`),colour="blue", lwd=1) +
  geom_ribbon(aes(ymin=simplesloper_q$`2.5%`, ymax=simplesloper_q$`97.5%`),fill="blue",alpha=0.2) +
  geom_line(aes(y=simpleslopep_q$`50%`),colour="red", lwd=1) +
  geom_ribbon(aes(ymin=simpleslopep_q$`2.5%`, ymax=simpleslopep_q$`97.5%`),fill="red",alpha=0.2) +
  labs(y = "Effect of SPI on TVP", x = "Z-Score of Percent Water Rights") + theme_bw(base_size = 11)
  

#Binomial Transformations

##return the antilogit of the intercept posterior marginal --> interpret as the average prob of y=1 (B&F) when all predictors at reference value or zero
Hint<-inla.tmarginal(function(x) exp(x)/(1+exp(x)),H$marginals.fixed[[1]])#build the transformed posterior
Hintquant<-inla.zmarginal(Hint) #obtain the summary stats for the intercept --> median prob of being B&F 

##return the exponentiated posterior of the coefficients --> interpret as the incremental change in the prob of B&F given a one unit increase in predictor
#H_probincr<-inla.emarginal(exp,H$marginals.fixed$pre_perc) #this is only for the mean, recommended to use the entire posterior to avoid biases
Hrip<-inla.tmarginal(function(x) exp(x),H$marginals.fixed$agrip_perc)#build the transformed posterior for riparian
Hripquant<-inla.zmarginal(Hrip)
Hrip_probincr<- (Hripquant$quant0.5-1) *100 #% change in the prob of B&F when Rip increases by 1 stdev
Hrip_probincr_025<- (Hripquant$quant0.025-1) *100 
Hrip_probincr_975<- (Hripquant$quant0.975-1) *100 

Hpre<-inla.tmarginal(function(x) exp(x),H$marginals.fixed$agpre_perc)
Hprequant<-inla.zmarginal(Hpre)
Hpre_probincr<- (Hprequant$quant0.5-1) *100 
Hpre_probincr_025<- (Hprequant$quant0.025-1) *100 
Hpre_probincr_975<- (Hprequant$quant0.975-1) *100 

Happ<-inla.tmarginal(function(x) exp(x),H$marginals.fixed$agapprop_perc)
Happquant<-inla.zmarginal(Happ)
Happ_probincr<- (Happquant$quant0.5-1) *100 
Happ_probincr_025<- (Happquant$quant0.025-1) *100 
Happ_probincr_975<- (Happquant$quant0.975-1) *100 

Hspi<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$spi)
Hspiquant<-inla.zmarginal(Hspi)
Hspi_probincr<- (Hspiquant$quant0.5-1) *100 
Hspi_probincr_025<- (Hspiquant$quant0.025-1) *100 
Hspi_probincr_975<- (Hspiquant$quant0.975-1) *100 

Hrs<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`agrip_perc:spi`)
Hrsquant<-inla.zmarginal(Hrs)
Hrs_probincr<- (Hrsquant$quant0.5-1) *100 
Hrs_probincr_025<- (Hrsquant$quant0.025-1) *100 
Hrs_probincr_975<- (Hrsquant$quant0.975-1) *100 
SPI_Rip1<-(Hspi_probincr/100)*(Hrs_probincr/100) *100 # % change in the prob of B&F when SPI increases by 1 stdev and Rip = 1
rsint<-inla.emarginal(exp,(H$marginals.fixed$spi+H$marginals.fixed$`agrip_perc:spi`))
rsint<-(1-rsint)*100

Hps<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`agpre_perc:spi`)
Hpsquant<-inla.zmarginal(Hps)
Hps_probincr<- (Hpsquant$quant0.5-1) *100 
Hps_probincr_025<- (Hpsquant$quant0.025-1) *100 
Hps_probincr_975<- (Hpsquant$quant0.975-1) *100 
SPI_Pre1<-(Hspi_probincr)/100*(Hps_probincr)/100  
psint<-inla.emarginal(exp,(H$marginals.fixed$spi+H$marginals.fixed$`agpre_perc:spi`))
psint<-(1-psint)*100

Has<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`agapprop_perc:spi`)
Hasquant<-inla.zmarginal(Has)
Has_probincr<- (Hasquant$quant0.5-1) *100 
Has_probincr_025<- (Hasquant$quant0.025-1) *100 
Has_probincr_975<- (Hasquant$quant0.975-1) *100 
SPI_App1<-(Hspi_probincr/100)*(Has_probincr/100) *100 
asint<-inla.emarginal(exp,(H$marginals.fixed$spi+H$marginals.fixed$`approp_perc:spi`))
asint<-(1-asint)*100

Ho<-inla.tmarginal(function(x) exp(x), H$marginals.fixed$`ag_perc`)
Hoquant<-inla.zmarginal(Ho)




