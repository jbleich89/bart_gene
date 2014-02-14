##run dynatree safely script:

setwd("C:/Users/jbleich/workspace/bart_gene/simulations_for_var_selection")
library(dynaTree)
source("dynatree_var_sel.R")

calc_prec_rec = function(true_vars, regression_vars){
  true_vars_found = intersect(true_vars, regression_vars)
  tps = length(true_vars_found)
  fps = length(setdiff(regression_vars, true_vars_found))
  fns = length(setdiff(true_vars, true_vars_found))
  if (tps + fps == 0){
    precision = 0
  } else {
    precision = tps / (tps + fps)
  }
  list(
    precision = precision,
    recall = tps / (tps + fns)
  )
}


num_replicates = 100
n = 250
ps = c(20, 100, 200, 500, 1000)
po_props = c(0.01, 0.05, 0.1, 0.2)
sigsqs = c(1, 5, 20)

##figure 1
p = 200 
po_prop = .01 
sigsq = 20

rep_results = array(NA, c(1, 2, num_replicates))

for (nr in 1 : num_replicates){
  cat("replicate #", nr, "\n")
  
  #generate linear model data
  X = matrix(rnorm(n * p), ncol = p)
  p0 = ceiling(p * po_prop)
  true_vars = 1 : p0
  beta_vec = c(rep(1, p0), rep(0, p - p0))
  error = rnorm(n, 0, sqrt(sigsq))
  y = as.numeric(X %*% beta_vec + error)
  X = as.data.frame(X)
  colnames(X) = seq(1, p)
  
  dynatree_vars = as.numeric(var_sel_dynaTree(as.matrix(X), y, n_particles = 3000))
  
  
  obj = calc_prec_rec(true_vars, dynatree_vars)
  rep_results[1, , nr] = c(obj$precision, obj$recall)
  print(rep_results[1, , nr])
  
  write.csv(rep_results[1, , nr], file = paste("sim_results/partial_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, "_nr_", nr, "_revised_dtree_only.csv", sep = ""))	
  
}

results = matrix(0, nrow = 1, ncol = 2)
rownames(results) = ("dynaTree") 
colnames(results) = c("precision", "recall")


#now dump results in
for (nr in 1 : num_replicates){
  results = results + rep_results[, , nr]
}
results = results / num_replicates

#calculate F1
F1s = 2 * results[, 1] * results[, 2] / (results[, 1] + results[, 2])
results = cbind(results, F1s)


#save results
write.csv(results, file = paste("../bart_gene/simulations_for_var_selection/sim_results/complete_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, "_revised_dtree_only.csv", sep = ""))

