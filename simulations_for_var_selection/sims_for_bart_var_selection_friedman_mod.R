library(MASS)
tryCatch(library(glmnet), error = function(e){install.packages("glmnet")}, finally = library(glmnet))


LAST_NAME = "kapelner"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

if (NOT_ON_GRID){
	setwd("C:/Users/Kapelner/workspace/CGMBART_GPL/")
} else {
	setwd("../CGMBART_GPL/")
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
num_replicates = 50
n = 250
ps = c(25, 100, 200, 500, 1000, 5000)
sigsqs = c(1, 10, 50, 100)

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

rep_results = array(NA, c(7, 2, num_replicates))

######replicate a few times
for (nr in 1 : num_replicates){
#generate friedman model data
	X = matrix(rnorm(n * p), ncol = p)
	X = as.data.frame(X)
	colnames(X) = seq(1, p)
	error = rnorm(n, 0, sqrt(sigsq))
	y = 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5] + error
	true_vars = 1 : 5
	
	
	#now build bart machine	
	bart_machine = build_bart_machine(X, y, num_trees = 1, num_burn_in = 2000, run_in_sample = FALSE)
	
	#do var selection with bart
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = ifelse(NOT_ON_GRID, TRUE, FALSE))
	bart_ptwise_vars = sort(as.numeric(bart_variables_select_obj$important_vars_pointwise))
	bart_simul_max_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_max))
	bart_simul_se_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_se))
	
	#do var selection with a CV-min-RMSE
	bart_cv_vars = var_selection_by_permute_response_cv(bart_machine)
	
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
	
	#do var selection with lasso
	lasso_matrix_vars = coef(cv.glmnet(as.matrix(X), y, alpha = 1 , nfolds = 5, type.measure = "mse"))
	lasso_matrix_vars = which(lasso_matrix_vars != 0) - 1
	lasso_matrix_vars = sort(lasso_matrix_vars[lasso_matrix_vars != 0]) #kill intercept if it exists
	
	
	#### now what did we get right?
	obj = calc_prec_rec(true_vars, bart_cv_vars)
	rep_results[1, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_ptwise_vars)
	rep_results[2, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars)
	rep_results[3, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_se_vars)
	rep_results[4, , nr] = c(obj$precision, obj$recall)
	if (p < n){
		obj = calc_prec_rec(true_vars, stepwise_backward_vars)
		rep_results[5, , nr] = c(obj$precision, obj$recall)
		obj = calc_prec_rec(true_vars, stepwise_forward_vars)
		rep_results[6, , nr] = c(obj$precision, obj$recall)		
	}
	obj = calc_prec_rec(true_vars, lasso_matrix_vars)
	rep_results[7, , nr ] = c(obj$precision, obj$recall)
	
}

results = matrix(0, nrow = 7, ncol = 2)
rownames(results) = c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "stepwise_backward",  "stepwise_forward", "lasso")
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
write.csv(results, file = paste("../bart_gene/simulations_for_var_selection/var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))	


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

