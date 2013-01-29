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

test_ols_methods_for_gene= function(gene_num){

  out = list()
  out[["gene"]] = colnames(gene_train)[gene_num]
  
  y = c(gene_train[ , gene_num], gene_cv[ , gene_num])
  X = rbind(tf_train, tf_cv)
  
  ## OLS& Step-wise OLS
  ##train on 90%

  train_data = data.frame(X, y)
  reg = lm(y ~ . , data = train_data)
  ols_predict = predict(object = reg, newdata = as.data.frame(tf_test))
  out[["ols_rmse"]] = sqrt(sum((ols_predict - gene_test[ , gene_num]) ^ 2) / length(ols_predict))
  
  step_reg = stepAIC(object = reg, direction = "backward" , trace = F)
  step_predict = predict(object = step_reg, newdata = as.data.frame(tf_test))
  out[["step_rmse"]] =   sqrt(sum((step_predict - gene_test[, gene_num]) ^ 2) / length(step_predict))
  out[["step_num_vars"]] = length(step_reg$coefficients) - 1 ##exclude one for intercept
                           
  ##Lasso
  lasso_fit = cv.glmnet(X, y, alpha = 1 , nfolds = 5, type.measure = "mse")
  names(lasso_fit)
  lasso_fit$nzero
  lasso_fit$glmnet.fit
  lasso_predict = predict(lasso_fit, newx = tf_test)
  
  out[["lasso_rmse"]] = sqrt(sum((lasso_predict - gene_test[, gene_num]) ^ 2) / length(lasso_predict))
  out[["lasso_num_vars"]] = lasso_fit$nzero[which(lasso_fit$lambda == lasso_fit$lambda.1se)] ##change to lambda.min
  out
}


mclapply(X = 1 : length(gene_names), FUN = test_ols_methods_for_gene, mc.cores = NUM_CORES)



