library(dplyr)
library(gstat)

dr <- "C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Kate_Emily\\Data\\"

data <- tbl_df(read.csv(paste(dr,"gw_timeser.csv",sep=""), stringsAsFactors=FALSE))

#182 observations, 2000.01 to 2014.12, monthly

quality <- apply(data, MARGIN=1, FUN = function(x) length(x[!is.na(x)]))
#drop wells with fewer than 14 obs?
#n = 14345 wells total
#n = 5989 with > 14 observations

#lat/lon