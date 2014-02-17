LAST_NAME = "jbleich"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

library(MASS)
# tryCatch(library(glmnet), error = function(e){install.packages("glmnet")}, finally = library(glmnet))
# tryCatch(library(randomForest), error = function(e){install.packages("glmnet")}, finally = library(randomForest))
tryCatch(library(dynaTree), error = function(e){install.packages("dynaTree")}, finally = library(dynaTree))
# tryCatch(library(bartMachine), error = function(e){install.packages("bartMachine")}, finally = library(bartMachine))
tryCatch(library(spikeslab), error = function(e){install.packages("spikeslab")}, finally = library(spikeslab))
#tryCatch(library(BayesTree), error = function(e){install.packages("BayesTree")}, finally = library(BayesTree))

options(error=recover)

##Loads
if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/Desktop/Dropbox/BART_gene/", sep = ""))
	library(bartMachine, lib.loc = .libPaths()[2])
} else {
	setwd("../shane_data")
# 	library(bartMachine, lib.loc = "~/R")
}

##set bart memory
# set_bart_machine_memory(3000)



# load("all_results.RData")
# load("all_validations.RData")
priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/workspace/bart_gene/analysis_with_tuning", sep = ""))
} else {
	setwd("../bart_gene/analysis_with_tuning")
}

source("dynatree_var_sel.R")
source("simulation_params.R")


# if (NOT_ON_GRID){
# 	setwd(paste("C:/Users/", LAST_NAME, "/workspace/CGMBART_GPL/", sep = ""))
# } else {
# 	setwd("../CGMBART_GPL")
# }

RF_ALPHA = 0.05


##OLS component
test_all_methods_for_gene = function(gene_num){
  out = list()
  rmse_mat = matrix(NA, nrow = 1, ncol = 13)
  num_vars_vec = matrix(NA, nrow = 1, ncol = 13)
  
  rownames(rmse_mat) = colnames(gene_train)[gene_num]
  rownames(num_vars_vec) = rownames(rmse_mat)
  simulation_names = c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso-CV", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full", "Rob-Best", "dynaTree", "spikeslab")
   colnames(rmse_mat) = simulation_names
   colnames(num_vars_vec) = simulation_names
  
  y = c(gene_train[ , gene_num], gene_cv[ , gene_num])
  X = rbind(tf_train, tf_cv)
  
  ## OLS& Step-wise OLS
  ##train on 90%
  y_test = gene_test[ , gene_num]
  
#   train_data = data.frame(X, y)
#   reg = lm(y ~ . , data = train_data)
#   ols_predict = predict(object = reg, newdata = as.data.frame(tf_test))
#   rmse_mat[ , "OLS"] = sqrt(sum((ols_predict - y_test) ^ 2) / length(y_test))
  
#   step_reg = stepAIC(object = reg, direction = "backward" , trace = F)
#   step_predict = predict(object = step_reg, newdata = as.data.frame(tf_test))
#   rmse_mat[ , "Stepwise"] = sqrt(sum((step_predict - y_test) ^ 2) / length(y_test))
#   num_vars_vec[ , "Stepwise"] = length(step_reg$coefficients) - 1 ##exclude one for intercept
  
  ##Lasso CV
#   lasso_cv_fit = cv.glmnet(X, y, alpha = 1 , nfolds = 5, type.measure = "mse")
#   lasso_cv_predict = predict(lasso_cv_fit, newx = tf_test)  
#   rmse_mat[ , "Lasso-CV"] = sqrt(sum((lasso_cv_predict - y_test) ^ 2) / length(y_test))
#   num_vars_vec[ , "Lasso-CV"] = lasso_cv_fit$nzero[which(lasso_cv_fit$lambda == lasso_cv_fit$lambda.1se)] ##change to lambda.min
  
  ##dynatree
  n_particles = 3000
  vars = var_sel_dynaTree(as.matrix(X), y, n_particles = n_particles)
  #dt_final = dynaTrees(X = as.matrix(X[,vars]), y = y, model = "linear", N = n_particles, verb = 1, XX = tf_test[,vars], R = 2)
  dt_final = dynaTree(X = as.matrix(X[,vars]), y = y, model = "linear", N = n_particles)
  #dt_preds = apply(dt_final$mean,1,mean)
  dt_preds = predict(dt_final, XX = tf_test[,vars], quants = F)$mean
  rmse_mat[ , "dynaTree"] = sqrt(sum((dt_preds - y_test) ^ 2) / length(y_test))
  num_vars_vec[ , "dynaTree"] = length(vars)
  
  ##spike slab
  spikeslab_mod = spikeslab(x = as.matrix(X), y = y, verbose = F) #90%
  spikeslab_preds = predict(spikeslab_mod, tf_test)
  rmse_mat[ , "spikeslab"] = sqrt(sum((spikeslab_preds$yhat.gnet - y_test) ^ 2) / length(y_test))
  num_vars_vec[ , "spikeslab"] = sum(spikeslab_mod$gnet > 0) 
  
  
  ##Need original training data now for Lasso, and BART, RF
#   y_train = gene_train[, gene_num]
#   y_cv = gene_cv[, gene_num]
#   y_test = gene_test[, gene_num]
   
  ##Train RF
# #   ntree = 500
#   p = ncol(tf_train)
#   rf = randomForest(x = tf_train , y = y_train, ntree = 500 ,importance = TRUE)
#   rf_zscore = importance(rf, type = 1 ,scale = TRUE)
#   rf_point_vars = which(rf_zscore > qnorm(1 - RF_ALPHA))
#   rf_simul_vars = which(rf_zscore > qnorm(1 - RF_ALPHA / p))
#   
#   ##CV RF
#   if (length(rf_point_vars) == 0){
#     point_err = sum((y_cv - mean(y_train))^2)
#   } else {
#     rf_pointwise = randomForest(x = as.matrix(tf_train[, rf_point_vars]), y = y_train, ntree = ntree)
#     y_hat = predict(rf_pointwise, newdata = tf_cv)
#     point_err = sum((y_cv - y_hat)^2)
#   }
#   if (length(rf_simul_vars) == 0){
#     simul_err = sum((y_cv - mean(y_train))^2)
#   } else {
#     rf_simul = randomForest(x = as.matrix(tf_train[, rf_simul_vars]), y = y_train, ntree = ntree)
#     y_hat = predict(rf_simul, newdata = tf_cv)
#     simul_err = sum((y_cv - y_hat)^2)
#   }
# 
# 
#   ##Get RMSE on test and number of vars-Test
#   if(simul_err < point_err){
# 	  fin_vars = rf_simul_vars
#   } else {
# 	  fin_vars = rf_point_vars 
#   }
#   if(length(fin_vars) == 0 ){
#     rmse_mat[,"RF"] = sqrt(sum((y_test-mean(y_train))^2)/length(y_test))
#     num_vars_vec[,"RF"] = 0
#   } else{
# 	if (simul_err < point_err){
# 		fin_rf_obj = rf_simul
# 	} else {
# 		fin_rf_obj = rf_pointwise
# 	}
#     yhat_rf = predict(fin_rf_obj, newdata = tf_test)
#     rmse_mat[,"RF"] = sqrt(sum((y_test-yhat_rf)^2)/length(y_test))
#     num_vars_vec[,"RF"] = length(fin_vars)
#   }
# 
#   
#   
#   ##Lasso
#   ##Train
#   lasso_train_fit = glmnet(tf_train, y_train, alpha = 1)
#   ##CV
#   lasso_predict_mat = predict(lasso_train_fit, newx = tf_cv)
#   L2_mat = (lasso_predict_mat - y_cv)^2
#   lasso_rmse_cv = sqrt(apply(L2_mat, 2, sum) / length(y_cv))
#   lambda = which.min(lasso_rmse_cv)
#   ##Test
#   lasso_predict = predict(lasso_train_fit, newx = tf_test)[, lambda]
#   rmse_mat[ , "Lasso"] = sqrt(sum((lasso_predict - y_test) ^ 2) / length(y_test))
#   num_vars_vec[ , "Lasso"] = lasso_train_fit$df[lambda] ##change to lambda.min
#   
#   ##NULL model
#   rmse_mat[ , "Null"] = sqrt(sum((y_test-mean(y_train))^2)/length(y_test))
#   num_vars_vec[ , "Null"] = 0
#   ##BART Runs
#   rmse_list = all_validations[[gene_num]]
#   validation_list = all_results[[gene_num]]
#   
#   ##Find best TFs##################
#   min_rmse = Inf
#   best_tfs = c()
#   for (c in cs){
#     for (method in METHODS){
#       if (rmse_list[[as.character(c)]][["20"]][["0.05"]][[method]] < min_rmse){
#         best_tfs = names(validation_list[[as.character(c)]][["20"]][["0.05"]][[method]])
#         min_rmse = rmse_list[[as.character(c)]][["20"]][["0.05"]][[method]]
#       }
#     }
#   }
#   ##Set data
#   
#   
#   if (length(best_tfs) > 0){
#     #build training data from the subset
#     X_train = as.data.frame(tf_train[, best_tfs])
# 	names(X_train) = seq(from = 1, to = ncol(X_train), by = 1)
#     
#     #now run BART model with 200 trees
#     bart_machine = build_bart_machine(X_train, y_train,
#       num_trees = NUM_TREES_FOR_EVAL, 
#       num_burn_in = NUM_BURN_IN, 
#       num_iterations_after_burn_in = NUM_ITER_AFTER, 
# 	  run_in_sample = FALSE,
#       verbose = FALSE)
#     
#     #predict on cv data set only with important tf's      			
#     test_data = as.data.frame(tf_test[, best_tfs])
# 	names(test_data) = seq(from = 1, to = ncol(test_data), by = 1)
# 	
#     predict_obj = bart_predict_for_test_data(bart_machine, test_data, y_test)
# 	rmse_mat[, "BART-Best"] = predict_obj$rmse
	
##DEFUNCT	
# 	#### Rob BART-best
# 	rbart = bart(x.train = X_train, 
# 		y.train = y_train, 
# 		x.test = test_data,
# 		ntree = NUM_TREES_FOR_EVAL, 
# 		nskip = 1000, 
# 		ndpost = 1000, 
# 		verbose = FALSE)
# 	yhat_rbart = rbart$yhat.test.mean
# 	rmse_mat[, "Rob-Best"] = sqrt(sum((y_test - yhat_rbart)^2) / length(yhat_rbart))
	
# 	#### OLS BART-Best
# 	mod = lm(y_train ~ ., X_train)
# 	y_hat = predict(mod, test_data)
# 	rmse_mat[, "OLS-BART-Best"] = sqrt(sum((y_test - y_hat)^2) / length(y_hat))
# 	
#     
#   } else {
#     L2_err = sum((y_test - mean(y_train))^2)
#     null_rmse = sqrt(L2_err / length(y_test))
# 	rmse_mat[, "BART-Best"] = null_rmse
# 	rmse_mat[, "Rob-Best"] = null_rmse
# 	rmse_mat[, "OLS-BART-Best"] = null_rmse
# 	
#   } 
#   
#   num_vars_vec[, "BART-Best"] = length(best_tfs)	
#   num_vars_vec[, "Rob-Best"] = length(best_tfs)
#   num_vars_vec[, "OLS-BART-Best"] = length(best_tfs)
#   
#   #####
#   
#   ##Find best TFs for simult. max#####
#   min_rmse = Inf
#   best_tfs_s_max = c()
#   for(c in cs){
#     if(rmse_list[[as.character(c)]][["20"]][["0.05"]][[METHODS[2]]] < min_rmse){
#       best_tfs_s_max = names(validation_list[[as.character(c)]][["20"]][["0.05"]][[METHODS[2]]])
#       min_rmse = rmse_list[[as.character(c)]][["20"]][["0.05"]][[METHODS[2]]]
#     }
#   }
#   
#   
#   if (length(best_tfs_s_max) > 0){
# 	#build training data from the subset
# 	X_train = as.data.frame(tf_train[, best_tfs_s_max])
#     
#     #now run BART model with 200 trees
#     bart_machine = build_bart_machine(X_train, y_train,
#       num_trees = NUM_TREES_FOR_EVAL, 
#       num_burn_in = NUM_BURN_IN, 
#       num_iterations_after_burn_in = NUM_ITER_AFTER, 
# 	  run_in_sample = FALSE,
#       verbose = FALSE)
#     
#     #predict on cv data set only with important tf's
# 	test_data = as.data.frame(tf_test[, best_tfs_s_max])
#     predict_obj = bart_predict_for_test_data(bart_machine, test_data, y_test)
# 	destroy_bart_machine(bart_machine)
#     bart_rmse_s_max = predict_obj$rmse						
#     
#   } else {
#     L2_err = sum((y_test - mean(y_train))^2)
#     bart_rmse_s_max = sqrt(L2_err / length(y_test))						
#   }
#   rmse_mat[ , "BART-S.Max"] = bart_rmse_s_max
#   num_vars_vec[, "BART-S.Max"] = length(best_tfs_s_max)
#   
#   #run a BART on the full dataset
# 
# 	bart_machine = build_bart_machine(as.data.frame(tf_train), y_train,
# 	  num_trees = NUM_TREES_FOR_EVAL, 
# 	  num_burn_in = NUM_BURN_IN, 
# 	  num_iterations_after_burn_in = NUM_ITER_AFTER, 
# 	  run_in_sample = FALSE,
# 	  verbose = FALSE)
# 
# 	predict_obj = bart_predict_for_test_data(bart_machine, data.frame(tf_test), y_test)
# 	bart_rmse = predict_obj$rmse
# 	rmse_mat[ , "BART-Full"] = bart_rmse
	  
	list(rmse_mat = rmse_mat, num_vars_vec = num_vars_vec)
}



 run_combined_tests = function(gene_num){
  all_methods = test_all_methods_for_gene(gene_num)
  rmse_results = all_methods$rmse_mat
  num_var_results = all_methods$num_vars_vec
  print(rmse_results)
  print(num_var_results)                        
  save(rmse_results, file = paste("rmse_results_dt_spike_", gene_num, ".RData", sep = ""))
  save(num_var_results, file = paste("num_var_results_dt_spike_", gene_num, ".RData", sep =""))
}
  

###set simulation parameters here
args = commandArgs(TRUE)
if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}
if (NOT_ON_GRID){
	iter_num = 997
	print(paste("iter_num:", iter_num))	
#	run_combined_tests(iter_num+1)
#	run_combined_tests(iter_num+2)
#	run_combined_tests(iter_num+3)
#	run_combined_tests(iter_num+4)
#	run_combined_tests(iter_num+5)
#	run_combined_tests(iter_num+6)
#	run_combined_tests(iter_num+7)
#	run_combined_tests(iter_num+8)
#	run_combined_tests(iter_num+9)
#	run_combined_tests(iter_num+10)
#	run_combined_tests(iter_num+11)
#	run_combined_tests(iter_num+12)	
}

print(iter_num)
gene_to_run = c(76,307,1230,1338,1552,1633,1779,1935,1991,2273,2671,2777,2939,3050,3440,3572,3705,3815,3853,3999,4124,4127,4130,4133,4135,4136,4137,4138,4139,4318,4338,4454,4822,4897,5013,5122,5444,5704,6026,2742,4040)

run_combined_tests(gene_to_run[iter_num])






 