library(MASS)
tryCatch(library(glmnet), error = function(e){install.packages("glmnet")}, finally = library(glmnet))
tryCatch(library(randomForest), error = function(e){install.packages("glmnet")}, finally = library(randomForest))



options(error = recover)

LAST_NAME = "kapelner"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

if (NOT_ON_GRID){
	setwd("C:/Users/Kapelner/workspace/bart_gene/simulations_for_var_selection")
} else {
	setwd("simulations_for_var_selection")
}

source("rf_cv_validator.R")

if (NOT_ON_GRID){
	setwd("C:/Users/Kapelner/workspace/CGMBART_GPL/")
} else {
	setwd("../../CGMBART_GPL/")
}

source("r_scripts/bart_package.R")
source("r_scripts/bart_package_variable_selection.R")


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

###
num_replicates = 25
n = 250
ps = c(25, 100, 200, 500, 1000)
sigsqs = c(1, 100, 625)

rf_alpha = .05

param_mat = as.data.frame(matrix(NA, nrow = length(ps) * length(sigsqs), ncol = 2))
colnames(param_mat) = c("p", "sigsq")
i = 1
for (p in ps){
	for (sigsq in sigsqs){
		param_mat[i, ] = c(p, sigsq)
		i = i + 1
	}
}


#read in arguments supplied by qsub - this will tell use which gene to process
args = commandArgs(TRUE)
print(paste("args:", args))

if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}
if (NOT_ON_GRID){
	iter_num = 1
	set_bart_machine_num_cores(4)
}

p = param_mat[iter_num, 1]
sigsq = param_mat[iter_num, 2]

rep_results = array(NA, c(18, 2, num_replicates))

######replicate a few times
for (nr in 1 : num_replicates){
	cat("replicate #", nr, "\n")
	#generate friedman model data
	X = matrix(rnorm(n * p), ncol = p)
	X = as.data.frame(X)
	colnames(X) = seq(1, p)
	error = rnorm(n, 0, sqrt(sigsq))
	y = 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5] + error
	true_vars = 1 : 5
	
	#do var selection with bart WITHOUT PRIOR
	bart_machine = build_bart_machine(X, y, num_trees = 1, num_burn_in = 2000, run_in_sample = FALSE, verbose = FALSE)	
	#now build bart machine WITHOUT PRIOR
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = ifelse(NOT_ON_GRID, TRUE, FALSE))
	bart_ptwise_vars = sort(as.numeric(bart_variables_select_obj$important_vars_pointwise))
	bart_simul_max_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_max))
	bart_simul_se_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_se))	
	#do var selection with a CV-min-RMSE
	bart_cv_vars = var_selection_by_permute_response_cv(bart_machine)$important_vars_cv
	destroy_bart_machine(bart_machine)
	
	#now build bart machine	WITH GOOD PRIOR
	bart_machine = build_bart_machine(X, y, 
			cov_prior_vec = c(rep(2, 5), rep(1, p - 5)), ### 5 is a magic number because there are always p0=5 variables that matter in the Friedman data
			num_trees = 1, 
			num_burn_in = 2000, 
			run_in_sample = FALSE, 
			verbose = FALSE)	
	#do var selection with bart WITH  PRIOR
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = ifelse(NOT_ON_GRID, TRUE, FALSE))
	bart_ptwise_vars_good_prior = sort(as.numeric(bart_variables_select_obj$important_vars_pointwise))
	bart_simul_max_vars_good_prior = sort(as.numeric(bart_variables_select_obj$important_vars_simul_max))
	bart_simul_se_vars_good_prior = sort(as.numeric(bart_variables_select_obj$important_vars_simul_se))	
	#do var selection with a CV-min-RMSE
	bart_cv_vars_good_prior = var_selection_by_permute_response_cv(bart_machine)$important_vars_cv
	destroy_bart_machine(bart_machine)
	
	
	#now build bart machine	WITH BAD PRIOR
	bart_machine = build_bart_machine(X, y, 
			cov_prior_vec = c(rep(1, p - 5), rep(2, 5)), ### 5 is a magic number because there are always p0=5 variables that matter in the Friedman data
			num_trees = 1, 
			num_burn_in = 2000, 
			run_in_sample = FALSE, 
			verbose = FALSE)	
	#do var selection with bart WITH PRIOR
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = ifelse(NOT_ON_GRID, TRUE, FALSE))
	bart_ptwise_vars_bad_prior = sort(as.numeric(bart_variables_select_obj$important_vars_pointwise))
	bart_simul_max_vars_bad_prior = sort(as.numeric(bart_variables_select_obj$important_vars_simul_max))
	bart_simul_se_vars_bad_prior = sort(as.numeric(bart_variables_select_obj$important_vars_simul_se))	
	#do var selection with a CV-min-RMSE
	bart_cv_vars_bad_prior = var_selection_by_permute_response_cv(bart_machine)$important_vars_cv
	destroy_bart_machine(bart_machine)
	
	#do var selection with stepwise
	if (p < n){
		step_reg = stepAIC(object = lm(y ~ ., data = X), direction = "backward" , trace = F)
		stepwise_backward_vars = names(step_reg$coefficients)
		stepwise_backward_vars = as.numeric(gsub("`", "", stepwise_backward_vars))
		stepwise_backward_vars = sort(stepwise_backward_vars[!is.na(stepwise_backward_vars)])
		
		step_reg = stepAIC(object = lm(y ~ ., data = X), direction = "forward" , trace = F)
		stepwise_forward_vars = names(step_reg$coefficients)
		stepwise_forward_vars = as.numeric(gsub("`", "", stepwise_forward_vars))
		stepwise_forward_vars = sort(stepwise_forward_vars[!is.na(stepwise_forward_vars)])		
	}
  
    ##do var selection with RF
    rf = randomForest(x = X , y = y, ntree = 500 ,importance = T)
    rf_zscore = importance(rf, type=1 ,scale=T)
	rf_point_vars = which(rf_zscore > qnorm(1 - rf_alpha))
	rf_simul_vars = which(rf_zscore > qnorm(1 - rf_alpha / p))
	rf_cv_vars = rf_cv_var_selection(X, y, 500, rf_alpha)$important_vars_cv

  
	#do var selection with lasso
	lasso_matrix_vars = coef(cv.glmnet(as.matrix(X), y, alpha = 1 , nfolds = 5, type.measure = "mse"))
	lasso_matrix_vars = which(lasso_matrix_vars != 0) - 1
	lasso_matrix_vars = sort(lasso_matrix_vars[lasso_matrix_vars != 0]) #kill intercept if it exists
	
	
	#### now what did we get right?
	#BART vanilla
	obj = calc_prec_rec(true_vars, bart_cv_vars)
	rep_results[1, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_ptwise_vars)
	rep_results[2, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars)
	rep_results[3, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_se_vars)
	rep_results[4, , nr] = c(obj$precision, obj$recall)
	#BART with good prior
	obj = calc_prec_rec(true_vars, bart_cv_vars_good_prior)
	rep_results[5, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_ptwise_vars_good_prior)
	rep_results[6, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars_good_prior)
	rep_results[7, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_se_vars_good_prior)
	rep_results[8, , nr] = c(obj$precision, obj$recall)	
	#BART with bad prior
	obj = calc_prec_rec(true_vars, bart_cv_vars_bad_prior)
	rep_results[9, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_ptwise_vars_bad_prior)
	rep_results[10, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars_bad_prior)
	rep_results[11, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_se_vars_bad_prior)
	rep_results[12, , nr] = c(obj$precision, obj$recall)	
	#the stepwises
	if (p < n){
		obj = calc_prec_rec(true_vars, stepwise_backward_vars)
		rep_results[13, , nr] = c(obj$precision, obj$recall)
		obj = calc_prec_rec(true_vars, stepwise_forward_vars)
		rep_results[14, , nr] = c(obj$precision, obj$recall)		
	}
	#LASSO
	obj = calc_prec_rec(true_vars, lasso_matrix_vars)
	rep_results[15, , nr] = c(obj$precision, obj$recall)
	#The RF methods
	obj = calc_prec_rec(true_vars, rf_cv_vars)
	rep_results[16, , nr] = c(obj$precision, obj$recall) 	
	obj = calc_prec_rec(true_vars, rf_point_vars)
	rep_results[17, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, rf_simul_vars)
	rep_results[18, , nr] = c(obj$precision, obj$recall)
 
		
}

results = matrix(0, nrow = 18, ncol = 2)
rownames(results) = c(
		"BART_CV", 
		"BART_pointwise", 
		"BART_simul_max", 
		"BART_simul_se", 
		"BART_CV_good_prior", 
		"BART_pointwise_good_prior", 
		"BART_simul_max_good_prior", 
		"BART_simul_se_good_prior", 
		"BART_CV_bad_prior", 
		"BART_pointwise_bad_prior", 
		"BART_simul_max_bad_prior", 
		"BART_simul_se_bad_prior", 		
		"stepwise_backward",  
		"stepwise_forward", 
		"lasso", 
		"RF_CV",
		"RF_point", 
		"RF_simul") 
colnames(results) = c("precision", "recall")

#now dump results in as averages
for (nr in 1 : num_replicates){
	results = results + rep_results[, , nr]
}
results = results / num_replicates

#calculate F1
F1s = 2 * results[, 1] * results[, 2] / (results[, 1] + results[, 2])
results = cbind(results, F1s)

#save results
write.csv(results, file = paste("../bart_gene/simulations_for_var_selection/sim_results/complete_var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))	


##load results and print them to xtable
#sink("all_results_friedman.tex")
#tryCatch(library(xtable), error = function(e){install.packages("xtable")}, finally = library(xtable))
#for (p in ps){
#	for (sigsq in sigsqs){
#		rep_results = read.csv(paste("var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))
#		print(xtable(rep_results, caption = paste("$p = ", p, ",~\\sigsq = ", sigsq, "$\n", sep ="")))
#	}
#}
#sink()

