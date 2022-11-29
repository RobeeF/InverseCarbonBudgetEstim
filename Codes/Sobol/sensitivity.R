#############
# Sensitivity
#############
root.path = "your_path_to_data/"# Modify with the path on your machine
setwd(root.path) # Modify with the path on your machine

source('Codes/Sobol/model.R')

# install.packages("TSP", dependencies = TRUE)
# install.packages('RANN')
# install.packages('sensitivity')

library(sensitivity)
library(ggplot2)
library(latex2exp)
library(cowplot)



anderson_labels = c(TeX('$\\psi$'), TeX('$\\omega_a$'), TeX('$\\omega_FL$'),
                    TeX('$\\alpha$'), TeX('$\\phi_v$'),
                    TeX('$\\beta_v$'), TeX('$K_v$'), TeX('$\\phi_v$'),
                    TeX('$\\beta_v$'), TeX('$K_v$'), TeX('$\\phi_z$'),
                    TeX('$\\beta_z$'), TeX('$\\lambda_z$'), TeX('$K_z$'),
                    TeX('$\\phi_h$'), TeX('$\\beta_H$'), TeX('$\\lambda_h$'),
                    TeX('$K_h$'), TeX('$\\zeta$'), TeX('$zi2$'))

str_labels = c('psi', 'w_A', 'w_FL', 'alpha',
               'vatt_DOC' ,'vatt_beta', 'vatt_npe',
               'vfl_DOC','vfl_beta', 'vfl_npe', 'z_DOC',
               'z_beta', 'z_D2', 'z_npe', 'h_DOC', 'h_beta',
               'h_D2', 'h_npe', 'zi', 'zi2'
)

#==================================================
# Sobol knn
#==================================================

N <- 20000
k <- 20 # without the CFs
X <- data.frame(matrix(runif(k * N), nrow = N))
# k2 <- 2 # The 2 CFs
# X2 = data.frame(matrix(runif(k2 * N, min = 0, max = 6) , nrow = N))
# X = cbind(X, X2)

stats_labels = c('Production_NonSinking', 'Production_Sinking', 'Respiration_Sinking', 'Respiration_zoo')


order.one.indices = matrix(0, length(stats_labels), k )
colnames(order.one.indices) <- str_labels

total.indices = matrix(0, length(stats_labels), k )
colnames(total.indices) <- str_labels


for (flow in 1:4){
  cat('flow number', flow, '\n')
  start = Sys.time()
  x <- sobolshap_knn(model = anderson, X = X, U = 1,
                     method = "rank", flux_out_nb = flow)

  colnames(x$S) <- str_labels
  order.one.indices[flow,] = x$S

  x2 <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
  tot.indices = tell(x2,x$y)

  colnames(x2$S) <- str_labels
  total.indices[flow,] = x2$S

  end = Sys.time()
  print(end-start)

}

rownames(order.one.indices) <- stats_labels
rownames(total.indices) <- stats_labels

write.csv(order.one.indices, 'Results/Sobol/1st_order_Sobol.csv')
write.csv(total.indices, 'Results/Sobol/total_Sobol.csv')

#===============================
# Plotting utility
#===============================


# First order
png('Results/Sobol/1st_order.png',  width = 750, height = 650)
par(mfrow = c(2,2))
for (flow in 1:4){
  plot(order.one.indices[flow,], main = paste('First order Sobol indices:', stats_labels[flow]),
       ylab = 'First order Sobol indices', xlab = '', las=2, xaxt='n',
       col="darkred")
  axis(1, at=1:(k + k2), labels = anderson_labels, las = 2)
}
dev.off()

# Total order
png('Results/Sobol/total.png',  width = 750, height = 650)
par(mfrow = c(2,2))
for (flow in 1:4){
  plot(total.indices[flow,], main = paste('Total Sobol indices:', stats_labels[flow]),
       ylab = 'Total Sobol indices', xlab = '', las=2, xaxt='n',
      col="darkred")
  axis(1, at=1:(k + k2), labels = anderson_labels, las = 2)
}
dev.off()
