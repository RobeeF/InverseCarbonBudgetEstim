#======================================
# Utilities definition
#======================================

# The weighted euclidian distance :
euc.dist <- function(x1, x2, w = 1) sqrt(rowSums((w *(x1 - x2)) ^ 2))

# The labels used for outputing
anderson_labels_txt = c('POC', 'DOC', 'psi', 'wA', 'wFL', 'alpha', 'CF_A', 
                        'CF_FL')

# Compute the normalized error between the predicted out flows and the actual outflows
pred.error <- function(X){
  models.params = X[1:8] # The 8 first parameters are the input parameters (see line 20 for the names)
  out.flows = as.numeric(X[9:11]) # The three last parameters are the values of the outflows
  pred = anderson(as.numeric(models.params))  
  w = 1 / as.numeric(out.flows) # Weight by the inverse of the flows to make all outflows comparable
  d = euc.dist(matrix(pred,1), matrix(out.flows,1), w)
  d
} # Watch out, the number of parameters is hard-coded for the moment in this function...
