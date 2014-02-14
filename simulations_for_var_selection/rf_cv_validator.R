
rf_cv_var_selection = function(X, y, ntree, rf_alpha, k_folds = 5, type = 1){
	n_cv = nrow(X)
	if (k_folds <= 1 || k_folds > n_cv){
		stop("The number of folds must be at least 2 and less than or equal to n, use \"Inf\" for leave one out")
	}
	
	
	if (k_folds == Inf){ #leave-one-out
		k_folds = n_cv
	}	
	
	holdout_size = round(n_cv / k_folds)
	split_points = seq(from = 1, to = n_cv, by = holdout_size)[1 : k_folds]
	
	L2_err_mat = matrix(NA, nrow = k_folds, ncol = 2)
	colnames(L2_err_mat) = c("pointwise", "simul")
	
	for (k in 1 : k_folds){
		cat("cv #", k, "\n", sep = "")
		#find out the indices of the holdout sample
		holdout_index_i = split_points[k]
		holdout_index_f = ifelse(k == k_folds, n_cv, split_points[k + 1] - 1)
		
		#pull out training and test data
		training_X_k = X[-c(holdout_index_i : holdout_index_f), ]
		training_y_k = y[-c(holdout_index_i : holdout_index_f)]
		test_X_k = X[holdout_index_i : holdout_index_f, ]
		text_y_k = y[holdout_index_i : holdout_index_f]
		
		#do the stuff
		rf = randomForest(x = as.matrix(training_X_k) , y = training_y_k, ntree = ntree, importance = T)
		rf_zscore = importance(rf, type = type ,scale = T)
		rf_point_vars = which(rf_zscore > qnorm(1 - rf_alpha))
		rf_simul_vars = which(rf_zscore > qnorm(1 - rf_alpha / p))
		
		#now get L2s
		if (length(rf_point_vars) == 0){
			L2_err_mat[k, "pointwise"] = sum((text_y_k - mean(training_y_k))^2)
		} else {
			rf_pointwise = randomForest(x = as.matrix(training_X_k[, rf_point_vars]), y = training_y_k, ntree = ntree)
			y_hat = predict(rf_pointwise, newdata = test_X_k)
			L2_err_mat[k, "pointwise"] = sum((text_y_k - y_hat)^2)
		}
		if (length(rf_simul_vars) == 0){
			L2_err_mat[k, "simul"] = sum((text_y_k - mean(training_y_k))^2)
		} else {
			rf_simul = randomForest(x = as.matrix(training_X_k[, rf_simul_vars]), y = training_y_k, ntree = ntree)
			y_hat = predict(rf_simul, newdata = test_X_k)
			L2_err_mat[k, "simul"] = sum((text_y_k - y_hat)^2)
		}
	}
	
	#now extract the lowest oos-L2 to find the "best" method for variable selection
	L2_err_by_method = colSums(L2_err_mat)
	min_var_selection_method = colnames(L2_err_mat)[which(L2_err_by_method == min(L2_err_by_method))]
	
	#now (finally) do var selection on the entire data and then return the vars from the best method found via cross-validation
	cat("final", "\n")
	rf = randomForest(x = X , y = y, ntree = ntree, importance = T)
	rf_zscore = importance(rf, type = type ,scale=T)
	rf_point_vars = which(rf_zscore > qnorm(1 - rf_alpha))
	rf_simul_vars = which(rf_zscore > qnorm(1 - rf_alpha / p))
	
	if (min_var_selection_method == "pointwise"){
		important_vars_cv = rf_point_vars
	} else {
		important_vars_cv = rf_simul_vars
	}
	
	list(best_method = min_var_selection_method, important_vars_cv = important_vars_cv)
}