#====================================================
#Barplot C budget
#====================================================


rm(list=ls())
cat("\014")


library(ggplot2)

# Paths to update on your machine:
# data.folder: Where to look for the outflows
# res.folder: Where to find the PGE/ CFs parameters estimates

root.folder = "your_path_to_data/" # Modify with the path on your machine
data.folder = paste0(root.folder, "Data/")
res.folder = paste0(root.folder, "Results/params_estimation/params")


flows = read.csv(file.path(data.folder, 'out_flows_4flows.csv'), sep = ';')
names(flows) = c('Cruise', 'Station', 'Depths',
                 'POC net', 'DOC net', 'prod non sinking',
                 'prod sinking', 'resp sinking', 'resp zoo')

Param = read.csv(file.path(res.folder, 'DY032 PAP config_CFs_0.5_PGE_Bounded_100runs .csv '), sep = ",")

# Get the sample location meta data
cruise = flows[1,'Cruise']
station = flows[1,'Station']
station = gsub('/', '-', station)
