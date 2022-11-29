# Integration using the trapeze method :

#***********************************
rm(list=ls())
cat("\014")
#**********************************


#working directory
root.folder = "your_path_to_data"
setwd(paste0(root.folder,"Data")

#Load data
df<-read.table("pour_power_law", header=T, sep="\t", na.strings = "NA", dec = ".")
head("df")


library(cubicBsplines)

trapeze(df$Depth_m, df$mgC_m_3d_1) # ! Be careful, the results is actually the cumsum.
