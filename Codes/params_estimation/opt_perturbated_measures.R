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

root.folder = "your_path_to_data" # Modify with the path on your machine
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

# From Anderson:
psi = 0.76
alpha = 0.5

#=========================================
# Stability of the optimisation process
#=========================================

nb.comb.params = 10
pge.bounds = 10
CF.bounds = 6
total.flow.nb = 6
cv.df = matrix(NA, total.flow.nb, nb_params - nb_fixed_entries)

for (i in 1:total.flow.nb){
  best.params.values = matrix(NA,nb.comb.params, nb_params - nb_fixed_entries)

  for (comb in 1:nb.comb.params){
    all.flows = as.numeric(c(poc, doc, c(out.flows[1,])))
    all.flows[i] = all.flows[i] * (1 + (runif(1) - 0.5) * 0.2)

    # The initial parameters values to launch the algorithm
    starting.values = c(all.flows[1], all.flows[2], runif(1),
                          runif(1)/pge.bounds, runif(1)/pge.bounds,
                          runif(1), 0.5, 0.5, all.flows[3],
                        all.flows[4], all.flows[5],
                        all.flows[6])
    starting.values = as.numeric(starting.values)

    # Lower/upper bounds contraints for each parameter.
    # The optimisation algorithm needs a non-null range for all params
    # => Use an epsilon to have a [poc - epsilon, poc + epsilon] range (same for the doc)
    epsilon = 1E-8

    lb = c(all.flows[1] - epsilon, all.flows[2] - epsilon,
           0, 0, 0, 0, 0.5 - epsilon, 0.5 - epsilon,
           all.flows[3] - epsilon,
           all.flows[4] - epsilon, all.flows[5] - epsilon,
           all.flows[6] - epsilon)
    lb = as.numeric(lb)
    ub = c(all.flows[1] + epsilon, all.flows[2] + epsilon,
           1, 1/pge.bounds, 1/pge.bounds, 1, 0.5 + epsilon,
           0.5 + epsilon, all.flows[3] + epsilon,
           all.flows[4] + epsilon, all.flows[5] + epsilon,
           all.flows[6] + epsilon)
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
    best.params.values[comb,] = best.comb
  }
  mean_ = apply(best.params.values, 2, mean)
  sd_ = apply(best.params.values, 2, sd)
  cv = sd_ / mean_
  cv.df[i,] = cv
}

res = data.frame(cv.df)
colnames(res) = labels
rownames(res) = c('poc', 'doc', 'prod non sinking', 'prod sinking',
                  'resp sinking', 'resp zoo')

res = round(res * 100, 2)
write.csv(res[, 1:4], file.path(res.folder, 'var_estim_var_fluxes.csv'))
