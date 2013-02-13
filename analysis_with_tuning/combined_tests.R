LAST_NAME = "kapelner"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

library(MASS)
tryCatch(library(glmnet), error = function(e){install.packages("glmnet")}, finally = library(glmnet))
tryCatch(library(randomForest), error = function(e){install.packages("randomForest")}, finally = library(randomForest))
tryCatch(library(BayesTree), error = function(e){install.packages("BayesTree")}, finally = library(BayesTree))

##need to make as argument
gene_num = 2

##Loads
if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/Desktop/Dropbox/BART_gene/", sep = ""))
} else {
	setwd("../shane_data")
}


load("all_results.RData")
load("all_validations.RData")
priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/workspace/bart_gene/analysis_with_tuning", sep = ""))
} else {
	setwd("../bart_gene/analysis_with_tuning")
}

source("simulation_params.R")


if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/workspace/CGMBART_GPL/", sep = ""))
} else {
	setwd("../CGMBART_GPL")
}

source("r_scripts/bart_package.R")


##OLS component
test_all_methods_for_gene = function(gene_num){
  out = list()
  rmse_mat = matrix(NA, nrow = 1, ncol = 11)
  num_vars_vec = matrix(NA, nrow = 1, ncol = 11)
  
  rownames(rmse_mat) = colnames(gene_train)[gene_num]
  rownames(num_vars_vec) = rownames(rmse_mat)
  simulation_names = c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso-CV", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full", "Rob-Best")
   colnames(rmse_mat) = simulation_names
   colnames(num_vars_vec) = simulation_names
  
  y = c(gene_train[ , gene_num], gene_cv[ , gene_num])
  X = rbind(tf_train, tf_cv)
  
  ## OLS& Step-wise OLS
  ##train on 90%
  y_test = gene_test[ , gene_num]
  
  train_data = data.frame(X, y)
  reg = lm(y ~ . , data = train_data)
  ols_predict = predict(object = reg, newdata = as.data.frame(tf_test))
  rmse_mat[ , "OLS"] = sqrt(sum((ols_predict - y_test) ^ 2) / length(y_test))
  
  step_reg = stepAIC(object = reg, direction = "backward" , trace = F)
  step_predict = predict(object = step_reg, newdata = as.data.frame(tf_test))
  rmse_mat[ , "Stepwise"] = sqrt(sum((step_predict - y_test) ^ 2) / length(y_test))
  num_vars_vec[ , "Stepwise"] = length(step_reg$coefficients) - 1 ##exclude one for intercept
  
  ##Lasso CV
  lasso_cv_fit = cv.glmnet(X, y, alpha = 1 , nfolds = 5, type.measure = "mse")
  lasso_cv_predict = predict(lasso_cv_fit, newx = tf_test)  
  rmse_mat[ , "Lasso-CV"] = sqrt(sum((lasso_cv_predict - y_test) ^ 2) / length(y_test))
  num_vars_vec[ , "Lasso-CV"] = lasso_cv_fit$nzero[which(lasso_cv_fit$lambda == lasso_cv_fit$lambda.1se)] ##change to lambda.min
  
  ##RF
  rf_fit = randomForest(x=X, y=y)
  rf_predict = predict(rf_fit , newdata = tf_test)
  rmse_mat[ , "RF"] = sqrt(sum((rf_predict - y_test) ^ 2) / length(y_test))
  
  
  ##Need original training data now for Lasso and BART
  y_train = gene_train[, gene_num]
  y_cv = gene_cv[, gene_num]
  y_test = gene_test[, gene_num]
  
  
  
  ##Lasso
  ##Train
  lasso_train_fit = glmnet(tf_train, y_train, alpha = 1)
  ##CV
  lasso_predict_mat = predict(lasso_train_fit, newx = tf_cv)
  L2_mat = (lasso_predict_mat - y_cv)^2
  lasso_rmse_cv = sqrt(apply(L2_mat, 2, sum) / length(y_cv))
  lambda = which.min(lasso_rmse_cv)
  ##Test
  lasso_predict = predict(lasso_train_fit, newx = tf_test)[, lambda]
  rmse_mat[ , "Lasso"] = sqrt(sum((lasso_predict - y_test) ^ 2) / length(y_test))
  num_vars_vec[ , "Lasso"] = lasso_train_fit$df[lambda] ##change to lambda.min
  
  ##NULL model
  rmse_mat[ , "Null"] = sqrt(sum((y_test-mean(y_train))^2)/length(y_test))
  num_vars_vec[ , "Null"] = 0
  ##BART Runs
  rmse_list = all_validations[[gene_num]]
  validation_list = all_results[[gene_num]]
  
  ##Find best TFs##################
  min_rmse = Inf
  best_tfs = c()
  for (c in cs){
    for (method in METHODS){
      if (rmse_list[[as.character(c)]][["20"]][["0.05"]][[method]] < min_rmse){
        best_tfs = names(validation_list[[as.character(c)]][["20"]][["0.05"]][[method]])
        min_rmse = rmse_list[[as.character(c)]][["20"]][["0.05"]][[method]]
      }
    }
  }
  ##Set data
  
  
  if (length(best_tfs) > 0){
    #build training data
    
    training_data = data.frame(tf_train, y = y_train)
    training_data = training_data[, c(best_tfs, "y")]
    
    #now run BART model with 200 trees
    bart_machine = build_bart_machine(training_data, 
      num_trees = NUM_TREES_FOR_EVAL, 
      num_burn_in = NUM_BURN_IN, 
      num_iterations_after_burn_in = NUM_ITER_AFTER, 
      verbose = FALSE)
    
    #predict on cv data set only with important tf's      			
    test_data = data.frame(tf_test, y = y_test)
    test_data = test_data[, c(best_tfs, "y")]
    predict_obj = bart_predict_for_test_data(bart_machine, test_data)
	rmse_mat[, "BART-Best"] = predict_obj$rmse
	
	
	#### Rob BART-best
	rbart = bart(x.train = training_data[, 1 : (ncol(training_data) - 1)], 
		y.train = y_train, 
		x.test = test_data[, 1 : (ncol(test_data) - 1)],
		ntree = NUM_TREES_FOR_EVAL, 
		nskip = 1000, 
		ndpost = 1000, 
		verbose = FALSE)
	yhat_rbart = rbart$yhat.test.mean
	rmse_mat[, "Rob-Best"] = sqrt(sum((test_data$y - yhat_rbart)^2) / length(yhat_rbart))
	
	#### OLS BART-Best
	mod = lm(y ~ ., training_data)
	y_hat = predict(mod, test_data)
	rmse_mat[, "OLS-BART-Best"] = sqrt(sum((test_data$y - y_hat)^2) / length(y_hat))
	
    
  } else {
    L2_err = sum((y_test - mean(y_train))^2)
    null_rmse = sqrt(L2_err / length(y_test))
	rmse_mat[, "BART-Best"] = null_rmse
	rmse_mat[, "Rob-Best"] = null_rmse
	rmse_mat[, "OLS-BART-Best"] = null_rmse
	
  } 
  
  num_vars_vec[, "BART-Best"] = length(best_tfs)	
  num_vars_vec[, "Rob-Best"] = length(best_tfs)
  num_vars_vec[, "OLS-BART-Best"] = length(best_tfs)
  
  #####
  
  ##Find best TFs for simult. max#####
  min_rmse = Inf
  best_tfs_s_max = c()
  for(c in cs){
    if(rmse_list[[as.character(c)]][["20"]][["0.05"]][[METHODS[2]]] < min_rmse){
      best_tfs_s_max = names(validation_list[[as.character(c)]][["20"]][["0.05"]][[METHODS[2]]])
      min_rmse = rmse_list[[as.character(c)]][["20"]][["0.05"]][[METHODS[2]]]
    }
  }
  
  
  if (length(best_tfs_s_max) > 0){
    #build training data
    
    training_data = data.frame(tf_train, y = y_train)
    training_data = training_data[, c(best_tfs_s_max, "y")]
    
    #now run BART model with 200 trees
    bart_machine = build_bart_machine(training_data, 
      num_trees = NUM_TREES_FOR_EVAL, 
      num_burn_in = NUM_BURN_IN, 
      num_iterations_after_burn_in = NUM_ITER_AFTER, 
      verbose = FALSE)
    
    #predict on cv data set only with important tf's        		
    test_data = data.frame(tf_test, y = y_test)
    test_data = test_data[, c(best_tfs_s_max, "y")]
    predict_obj = bart_predict_for_test_data(bart_machine, test_data)
	destroy_bart_machine(bart_machine)
    bart_rmse_s_max = predict_obj$rmse						
    
  } else {
    L2_err = sum((y_test - mean(y_train))^2)
    bart_rmse_s_max = sqrt(L2_err / length(y_test))						
  }
  rmse_mat[ , "BART-S.Max"] = bart_rmse_s_max
  num_vars_vec[, "BART-S.Max"] = length(best_tfs_s_max)
  
  #run a BART on the full dataset

  	training_data = data.frame(tf_train, y = y_train)
  	bart_machine = build_bart_machine(training_data, 
	  num_trees = NUM_TREES_FOR_EVAL, 
	  num_burn_in = NUM_BURN_IN, 
	  num_iterations_after_burn_in = NUM_ITER_AFTER, 
	  verbose = FALSE)
	test_data = data.frame(tf_test, y = y_test)
	predict_obj = bart_predict_for_test_data(bart_machine, test_data)
	bart_rmse = predict_obj$rmse
	rmse_mat[ , "BART-Full"] = bart_rmse
  
	list(rmse_mat = rmse_mat, num_vars_vec = num_vars_vec)
}



run_combined_tests = function(gene_num){
  all_methods = test_all_methods_for_gene(gene_num)
  rmse_results = all_methods$rmse_mat
  num_var_results = all_methods$num_vars_vec
  print(rmse_results)
  print(num_var_results)                        
  save(rmse_results, file = paste("rmse_results_", gene_num, ".RData", sep = ""))
  save(num_var_results, file = paste("num_var_results_", gene_num, ".RData", sep =""))
}
  

###set simulation parameters here
args = commandArgs(TRUE)
if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}
if (NOT_ON_GRID){
	iter_num = 25
	print(paste("iter_num:", iter_num))	
	run_combined_tests(iter_num+1)
	run_combined_tests(iter_num+2)
	run_combined_tests(iter_num+3)
	run_combined_tests(iter_num+4)
	run_combined_tests(iter_num+5)
	run_combined_tests(iter_num+6)
	run_combined_tests(iter_num+7)
	run_combined_tests(iter_num+8)
	run_combined_tests(iter_num+9)
	run_combined_tests(iter_num+10)
	run_combined_tests(iter_num+11)
	run_combined_tests(iter_num+12)	
}

run_combined_tests(iter_num)






 