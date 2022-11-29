#====================================================
# Optimisation : Anderson Model - 4 outflows version
#====================================================


rm(list=ls())
cat("\014")

#install.packages('dfoptim')
library(dfoptim)
library(stringr)
library(ggplot2)
library(cowplot)
library(tidyr)

# Paths to update on your machine:
# data.folder: Where to look for the outflows to which model outputs will be compared
# res.folder: Where to write the PGE/ CFs parameters estimates
# code.folder: The path to the Anderson model

root.folder = "your_path_to_data/" # Modify with the path on your machine
data.folder = paste0(root.folder, "Data/")
res.folder = paste0(root.folder, "Results/params_estimation")
code.folder = paste0(root.folder, "Codes/params_estimation")

source(file.path(code.folder, 'model.R'))
source(file.path(code.folder, 'utilities.R'))

#=========================================
# Data handling
#=========================================

nb_params = length(anderson_labels_txt) # The number of parameters to estimate
nb_fixed_entries = length(c('poc', 'doc')) # The in situ poc and doc measures
labels = anderson_labels_txt[(nb_fixed_entries + 1):nb_params]

# Fetch the observed outflows
flows = read.csv(file.path(data.folder, 'out_flows_4flows.csv'), sep = ';')
names(flows) = c('Cruise', 'Station', 'Depths',
                     'POC net', 'DOC net', 'prod non sinking',
                     'prod sinking', 'resp sinking', 'resp zoo')

out.flows = flows[,c('prod non sinking', 'prod sinking',
                     'resp sinking', 'resp zoo')]
nb_stats = length(out.flows)

# Get the poc and doc flows
poc = flows[1,'POC net']
doc = flows[1,'DOC net']

# Get the sample location meta data
cruise = flows[1,'Cruise']
station = flows[1,'Station']
station = gsub('/', '-', station)


#=========================================
# Stability of the optimisation process
#=========================================

nb.comb.params = 10
best.params.values = matrix(NA,nb.comb.params, nb_params - nb_fixed_entries + 1)
pge.bounds = 10
CF.bounds = 6

for (comb in 1:nb.comb.params){
  # The initial parameters values to launch the algorithm
  starting.values = c(c(poc, doc, runif(1), runif(1)/pge.bounds, runif(1)/pge.bounds, runif(1), 0.5, 0.5), out.flows[1,])
  #starting.values = c(c(poc, doc, 0.75, 0.02, 0.08, 0.5, 0.44, 0.44), out.flows[1,])
  starting.values = as.numeric(starting.values)

  # Lower/upper bounds contraints for each parameter.
  # The optimisation algorithm needs a non-null range for all params
  # => Use an epsilon to have a [poc - epsilon, poc + epsilon] range (same for the doc)
  epsilon = 1E-8

  lb = c(c(poc - epsilon, doc - epsilon, 0, 0, 0, 0, 0.5 - epsilon, 0.5 - epsilon), as.numeric(out.flows[1,]) - epsilon)
  lb = as.numeric(lb)
  ub = c(c(poc + epsilon, doc + epsilon, 1, 1/pge.bounds, 1/pge.bounds, 1, 0.5 + epsilon, 0.5 + epsilon), as.numeric(out.flows[1,]) + epsilon) # Lower bounds contraints for each parameter
  ub = as.numeric(ub)

  # Performs a gradient-free optimisation
  # The goal is to minimize the distance existing between the measured
  # outflows and their model counterparts
  start = Sys.time()
  res.station = nmkb(par = starting.values, fn = pred.error, lower=lb, upper=ub)
  end  = Sys.time()
  print(end - start)
  print(res.station$message)

  # Store the values of the parameters and the associated loss
  best.comb = res.station$par[(nb_fixed_entries + 1):nb_params]
  #best.comb = res.station$par[(nb_fixed_entries + 1):nb_params]

  best.comb = c(best.comb, res.station$value)
  best.params.values[comb,] = best.comb
}

res = data.frame(best.params.values)
colnames(res) = c(labels, 'loss')

#===========================================
# Find the best combinaison and compute CIs
#===========================================

# Estimation
#min.loss = min(res$loss)
#best.comb = res[res$loss == min.loss,]

best.comb = apply(res, 2, quantile, 0.5)
# CI bounds
q05.comb = apply(res, 2, quantile, 0.05) # Lower bound
q95.comb = apply(res, 2, quantile, 0.95) # Upper bound

# Store the results
estim = as.data.frame(t(rbind(best.comb, q05.comb, q95.comb)))
colnames(estim) = c('estimation', 'lb', 'ub')
estim = round(estim, 3)
#estim$estimation = formatC(estim$estimation, format = "e", digits = 1)
#estim$lb = formatC(estim$lb, format = "e", digits = 1)
#estim$ub = formatC(estim$ub, format = "e", digits = 1)

estim$CI = str_c('(', str_c(estim[,c('lb')],
                                        ' ; ', estim[, c('ub')] ), ')')
estim = subset(estim, select = -c(lb, ub) )

# Flows
models.params = c(c(poc, doc), estim[1:6,1], as.numeric(out.flows[1,])) # The 8 first parameters are the input parameters (see line 20 for the names)
pred = anderson(as.numeric(models.params))
(pred - out.flows) / out.flows

#===========================================
# Store and plot the results
#===========================================
estim.csv <- tibble::rownames_to_column(estim, "variable")

write.table(estim.csv,
            col.names = T,
            file.path(res.folder,
                      paste(cruise, station, 'config_PGE_bounded_CF_bounded',".csv")),
            row.names = F, sep = ',')
