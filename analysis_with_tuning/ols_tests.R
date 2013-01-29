##Test Script
library(MASS)
library(glmnet)
library(parallel)


setwd("~/Research_Genomics/")

priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

source("~/Research_Genomics/bart_gene/analysis_with_tuning/simulation_params_JB.R")

NUM_CORES=20

test_ols_methods_for_gene = function(gene_num){
  out = list()
  rmse_mat = matrix(NA,nrow=1,ncol=3)
  num_var_mat = matrix(NA,nrow=1,ncol=3)
  
  rownames(rmse_mat) = colnames(gene_train)[gene_num]
  rownames(num_var_mat) = rownames(rmse_mat)
  colnames(rmse_mat) = c("OLS","Stepwise","Lasso")
  colnames(num_var_mat) = colnames(rmse_mat)
  
  y = c(gene_train[ , gene_num], gene_cv[ , gene_num])
  X = rbind(tf_train, tf_cv)
  
  ## OLS& Step-wise OLS
  ##train on 90%
  y_test = gene_test[ , gene_num]
  
  train_data = data.frame(X, y)
  reg = lm(y ~ . , data = train_data)
  ols_predict = predict(object = reg, newdata = as.data.frame(tf_test))
  rmse_mat[ , "OLS"] = sqrt(sum((ols_predict - y_test) ^ 2) / length(ols_predict))
  num_var_mat[ , "OLS"] = 39
  
  step_reg = stepAIC(object = reg, direction = "backward" , trace = F)
  step_predict = predict(object = step_reg, newdata = as.data.frame(tf_test))
  rmse_mat[ , "Stepwise"] = sqrt(sum((step_predict - y_test) ^ 2) / length(step_predict))
  num_var_mat[ , "Stepwise"] = length(step_reg$coefficients) - 1 ##exclude one for intercept
  
  ##Lasso
  lasso_fit = cv.glmnet(X, y, alpha = 1 , nfolds = 5, type.measure = "mse")
  lasso_predict = predict(lasso_fit, newx = tf_test)
  
  rmse_mat[ , "Lasso"] = sqrt(sum((lasso_predict - y_test) ^ 2) / length(lasso_predict))
  num_var_mat[ , "Lasso"] = lasso_fit$nzero[which(lasso_fit$lambda == lasso_fit$lambda.1se)] ##change to lambda.min
  list("rmse" = rmse_mat, "num_vars" = num_var_mat)
}



mclapply(X = 1, FUN = test_ols_methods_for_gene, mc.cores = NUM_CORES)



