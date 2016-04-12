library(FactoMineR)
library(missMDA)
library(factoextra)
library(corrplot)
library(dplyr)


fn = "D:\\mca_data.csv"
ds <- tbl_df(read.csv(fn, stringsAsFactors = FALSE))

#for now predicting HAVE practiced
df <- data.frame(tvp = ds$px_tvp_y5, lulc = ds$lu07)
df <- data.frame(lapply(df, as.factor), stringsAsFactors = T)
df <- df[df$tvp != ]

#since tvp non categorical, won't work

#nb <- estim_ncpMCA(as.data.frame(ADP1), ncp.max = 5)
res.mca <- MCA(as.data.frame(df), graph = T, na.method = "NA")

#data display
#http://www.sthda.com/english/wiki/multiple-correspondence-analysis-essentials-interpretation-and-application-to-investigate-the-associations-between-categories-of-multiple-qualitative-variables-r-software-and-data-mining
dim <- dimdesc(res.mca)
summary(res.mca, nb.dec = 2, ncp = 3, nbelements = 25)  #file = my_file.txt to export
ev <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca)
var<- get_mca_var(res.mca)

plot(res.mca, choix = "var")
#corr plot
corrplot(var$contrib, is.corr = F)

#contribution of variables to dim 1
fviz_contrib(res.mca, choice = "var", axes = 1, top = 20) +   
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20), axis.text = element_text(size = 20), 
                     axis.title = element_text(size = 20), plot.title = element_text(size = 20))
               

#individual contribution
fviz_mca_ind(res.mca, label = "none", habillage = ADP1$major, addEllipses = T, ellipse.level = 0.95)

#ellipse
plotellipses(res.mca, keepvar = 1:4)
