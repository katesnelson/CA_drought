#################################################
#KRIG DIAGNOSTICS
##################################################

stplot(pred3)

kd<-as(pred3, "data.frame") #transform the stfdf to a data frame

pred_time_trend <- kd%>% group_by(time) %>% summarize(mean(var1.pred)) 
plot(x=pred_time_trend$time, y=pred_time_trend$`mean(var1.pred)`) #plot mean of predicted elev across space over time

pred_space_trend <- kd%>% group_by(sp.ID) %>% summarize(mean(var1.pred)) 
ksp<- left_join(kd, pred_space_trend, by = "sp.ID")
coordinates(ksp)<- ~x1-x2
ksp@proj4string <- CRS('+init=epsg:3310')
gridded(ksp)=TRUE
plot(ksp[ ,"mean(var1.pred)"]) #map mean of predicted elev across time over space

####################################################
#CHECK RMSE RESULTS for n =1000
####################################################

pred3_rmse$elev_obs <- as.numeric(as.character(pred3_rmse$elev_obs)) #convert from factor to numeric
pred3_rmse$elev_pred <- as.numeric(as.character(pred3_rmse$elev_pred))#convert from factor to numeric
pred3_rmse <- mutate(pred3_rmse, elev_diff = elev_obs - elev_pred)  #add a column for diff between observed and predicted

rmse3

rmse_time_trend <- pred3_rmse%>% group_by(V5) %>% summarize(mean(elev_diff)) 
plot(rmse_time_trend) #plot mean of predicted elev across space over time

rmse_space_trend <- pred3_rmse %>% group_by(lon, lat) %>% summarize(mean(elev_diff)) 
rmse_space_trend$lon <- as.numeric(as.character(rmse_space_trend$lon))
rmse_space_trend$lat <- as.numeric(as.character(rmse_space_trend$lat ))
rmse_space_trend<-as.data.frame(rmse_space_trend)
rmse_space_trend <- dplyr::rename(rmse_space_trend,  mean_elev_diff = `mean(elev_diff)`)

coordinates(rmse_space_trend)<- ~lon+lat
rmse_space_trend@proj4string <- CRS('+init=epsg:3310')
spplot(rmse_space_trend, main="mean_elev_diff")#map mean of elev residuals across time over space

##########################################
###RMSE RESULTS for nearly full HoldOUt###
##########################################

res$elev_obs <- as.numeric(as.character(res$elev_obs)) #convert from factor to numeric
res$elev_pred <- as.numeric(as.character(res$elev_pred))#convert from factor to numeric
res <- mutate(res, elev_diff = elev_obs - elev_pred)  #add a column for diff between observed and predicted
res$time_pred <- as.Date(as.character(res$V5))

rmse<-rmse(res$elev_obs, res$elev_pred)
rmse

res_time_trend <- res%>% group_by(V5) %>% summarize(mean(elev_diff)) 
plot(res_time_trend) #plot mean of elev residuals across space over time

res_space_trend <- res %>% group_by(lon, lat) %>% summarize(mean(elev_diff)) 
res_space_trend$lon <- as.numeric(as.character(res_space_trend$lon))
res_space_trend$lat <- as.numeric(as.character(res_space_trend$lat ))
res_space_trend<-as.data.frame(res_space_trend)
res_space_trend <- dplyr::rename(res_space_trend,  mean_elev_diff = `mean(elev_diff)`)

coordinates(res_space_trend)<- ~lon+lat
res_space_trend@proj4string <- CRS('+init=epsg:3310')
spplot(res_space_trend, main="mean_elev_diff")#map mean of  elev residuals across time over space




#plot RMSE value at a point in space (held out) through all time obs
time_obs <- held_out_sample%>% group_by(TIME) %>% summarize(mean(ELEV))
time_orig <-dfs_clean %>%group_by(TIME) %>% summarize(mean(ELEV))
time_pred <- res%>% group_by(time_pred) %>% summarize(mean(elev_pred))
time_m <- res %>% group_by(time_pred) %>% summarize(mean(elev_diff))
time_sd <- res %>% group_by(time_pred) %>% summarize(sd(elev_diff))
time_count <- held_out_sample%>% group_by(TIME) %>% summarize(count=n()/100)

plot(time_m, type="l", main="obs-pred", col="black", ylim=c(-10, 10))
lines((time_count), type='l', col='green')
lines(time_orig, type = 'l', lty = 2)
lines(time_obs, type='l', col='blue', lty=2)
lines(time_pred, type='l', col='red', lty=2)
legend("topright", c("obs", "pred"), lty=c(1,1), lwd=c(2.5,2.5), col=c("blue","red")) 


