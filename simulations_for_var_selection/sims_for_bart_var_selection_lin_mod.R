library(MASS)
tryCatch(library(glmnet), error = function(e){install.packages("glmnet")}, finally = library(glmnet))
tryCatch(library(randomForest), error = function(e){install.packages("glmnet")}, finally = library(randomForest))
tryCatch(library(dynaTree), error = function(e){install.packages("dynaTree")}, finally = library(dynaTree))
tryCatch(library(spikeslab), error = function(e){install.packages("spikeslab")}, finally = library(spikeslab))

options(error = recover)

LAST_NAME = "kapelner"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

if (NOT_ON_GRID){
	setwd("C:/Users/Kapelner/workspace/bart_gene/simulations_for_var_selection")
	library(bartMachine, lib.loc = .libPaths()[2])
} else {
	setwd("simulations_for_var_selection")
	library(bartMachine, lib.loc = "~/R")
}

source("rf_cv_validator.R")
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

###
NUM_ALGOS = 20

num_replicates = 50
n = 250
ps = c(20, 100, 200, 500, 1000)
po_props = c(0.01, 0.05, 0.1, 0.2)
sigsqs = c(1, 5, 20)

rf_alpha = .05

param_mat = as.data.frame(matrix(NA, nrow = length(ps) * length(po_props) * length(sigsqs), ncol = 3))
colnames(param_mat) = c("p", "po_prop", "sigsq")

i = 1
for (p in ps){
	for (po_prop in po_props){
		for (sigsq in sigsqs){
			param_mat[i, ] = c(p, po_prop, sigsq)
			i = i + 1
		}			
	}
}

param_mat = param_mat[c(26, 27, 32, 33, 40, 41, 43, 44), ]

#read in arguments supplied by qsub - this will tell use which gene to process
args = commandArgs(TRUE)
print(paste("args:", args))

if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}
if (NOT_ON_GRID){
	iter_num = 2
	set_bart_machine_num_cores(4)
}
##set bart memeory
set_bart_machine_memory(4000)

p = param_mat[iter_num, 1]
po_prop = param_mat[iter_num, 2]
sigsq = param_mat[iter_num, 3]


rep_results = array(NA, c(NUM_ALGOS, 2, num_replicates))




######replicate a few times
for (nr in num_replicates : 1){
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
	
	
	#now build bart machine WITHOUT PRIOR	  
	bart_machine = build_bart_machine(X, y, num_trees = 1, num_burn_in = 2000, run_in_sample = FALSE, verbose = FALSE)
	#do var selection with bart
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = FALSE)
# 	bart_ptwise_vars = sort(as.numeric(bart_variables_select_obj$important_vars_local_names))
	bart_simul_max_vars = sort(as.numeric(bart_variables_select_obj$important_vars_global_max_names))
# 	bart_simul_se_vars = sort(as.numeric(bart_variables_select_obj$important_vars_global_se_names))	
	#do var selection with a CV-min-RMSE
# 	bart_cv_vars = as.numeric(var_selection_by_permute_response_cv(bart_machine)$important_vars_cv)
	destroy_bart_machine(bart_machine)
	
	
	#now build bart machine WITH GOOD PRIOR	
	bart_machine = build_bart_machine(X, y, 
			cov_prior_vec = c(rep(2, p0), rep(1, p - p0)),
			num_trees = 1, 
			num_burn_in = 2000, 
			run_in_sample = FALSE, 
			verbose = FALSE)
	#do var selection with bart
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = FALSE)
# 	bart_ptwise_vars_good_prior = sort(as.numeric(bart_variables_select_obj$important_vars_local_names))
	bart_simul_max_vars_good_prior = sort(as.numeric(bart_variables_select_obj$important_vars_global_max_names))
# 	bart_simul_se_vars_good_prior = sort(as.numeric(bart_variables_select_obj$important_vars_global_se_names))	
	#do var selection with a CV-min-RMSE
# 	bart_cv_vars_good_prior = as.numeric(var_selection_by_permute_response_cv(bart_machine)$important_vars_cv)
 	destroy_bart_machine(bart_machine)
	
	#now build bart machine WITH BAD PRIOR	
	bart_machine = build_bart_machine(X, y, 
			cov_prior_vec = c(rep(1, p - p0), rep(2, p0)),
			num_trees = 1, 
			num_burn_in = 2000, 
			run_in_sample = FALSE, 
			verbose = FALSE)
	#do var selection with bart
	bart_variables_select_obj = var_selection_by_permute_response_three_methods(bart_machine, plot = FALSE)
# 	bart_ptwise_vars_bad_prior = sort(as.numeric(bart_variables_select_obj$important_vars_local_names))
	bart_simul_max_vars_bad_prior = sort(as.numeric(bart_variables_select_obj$important_vars_global_max_names))
# 	bart_simul_se_vars_bad_prior = sort(as.numeric(bart_variables_select_obj$important_vars_global_se_names))	

	#do var selection with a CV-min-RMSE
# 	bart_cv_vars_bad_prior = as.numeric(var_selection_by_permute_response_cv(bart_machine)$important_vars_cv)
	destroy_bart_machine(bart_machine)	
	
# 	#do var selection with stepwise
# 	if (p < n){
# 		step_reg = stepAIC(object = lm(y ~ ., data = X), direction = "backward" , trace = F)
# 		stepwise_backward_vars = names(step_reg$coefficients)
# 		stepwise_backward_vars = as.numeric(gsub("`", "", stepwise_backward_vars))
# 		stepwise_backward_vars = sort(stepwise_backward_vars[!is.na(stepwise_backward_vars)])
# 		
# 		step_reg = stepAIC(object = lm(y ~ ., data = X), direction = "forward" , trace = F)
# 		stepwise_forward_vars = names(step_reg$coefficients)
# 		stepwise_forward_vars = as.numeric(gsub("`", "", stepwise_forward_vars))
# 		stepwise_forward_vars = sort(stepwise_forward_vars[!is.na(stepwise_forward_vars)])		
# 	}
#   
# 	##do var selection with RF
# 	rf = randomForest(x = X , y = y, ntree = 500 ,importance = T)
# 	rf_zscore = importance(rf, type=1 ,scale=T)
# 	rf_point_vars = which(rf_zscore > qnorm(1 - rf_alpha))
# 	rf_simul_vars = which(rf_zscore > qnorm(1 - rf_alpha / p))  
# 	rf_cv_vars = rf_cv_var_selection(X, y, 500, rf_alpha)$important_vars_cv
# 	
# 	
# 	#do var selection with lasso
# 	lasso_matrix_vars = coef(cv.glmnet(as.matrix(X), y, alpha = 1 , nfolds = 5, type.measure = "mse"))
# 	lasso_matrix_vars = which(lasso_matrix_vars != 0) - 1
# 	lasso_matrix_vars = sort(lasso_matrix_vars[lasso_matrix_vars != 0]) #kill intercept if it exists
#   
#   
#   ##do var selection with dynaTree
#   dynatree_vars = as.numeric(var_sel_dynaTree(as.matrix(X), y, n_particles = 5000))
#   
#   ##do var selection with spikeslab
# 	spikeslab_mod = spikeslab(x = as.matrix(X), y = y, verbose = F, 
#                              bigp.smalln = ifelse(ncol(X) > nrow(X), T, F), max.var = 100)
#   spikeslab_vars = which(spikeslab_mod$gnet > 0)
  
	#### now what did we get right?
	#BART vanilla
# 	obj = calc_prec_rec(true_vars, bart_cv_vars)
# 	rep_results[1, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, bart_ptwise_vars)
# 	rep_results[2, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars)
	rep_results[3, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, bart_simul_se_vars)
# 	rep_results[4, , nr] = c(obj$precision, obj$recall)
# 	#BART with good prior
# 	obj = calc_prec_rec(true_vars, bart_cv_vars_good_prior)
# 	rep_results[5, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, bart_ptwise_vars_good_prior)
# 	rep_results[6, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars_good_prior)
 	rep_results[7, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, bart_simul_se_vars_good_prior)
# 	rep_results[8, , nr] = c(obj$precision, obj$recall)	
# 	#BART with BAD prior
# 	obj = calc_prec_rec(true_vars, bart_cv_vars_bad_prior)
# 	rep_results[9, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, bart_ptwise_vars_bad_prior)
# 	rep_results[10, , nr] = c(obj$precision, obj$recall)
 	obj = calc_prec_rec(true_vars, bart_simul_max_vars_bad_prior)
 	rep_results[11, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, bart_simul_se_vars_bad_prior)
# 	rep_results[12, , nr] = c(obj$precision, obj$recall)		
# 	#the stepwises
# 	if (p < n){
# 		obj = calc_prec_rec(true_vars, stepwise_backward_vars)
# 		rep_results[13, , nr] = c(obj$precision, obj$recall)
# 		obj = calc_prec_rec(true_vars, stepwise_forward_vars)
# 		rep_results[14, , nr] = c(obj$precision, obj$recall)		
# 	}
# 	#LASSO
# 	obj = calc_prec_rec(true_vars, lasso_matrix_vars)
# 	rep_results[15, , nr] = c(obj$precision, obj$recall)
#   
# 	#The RF methods
# 	obj = calc_prec_rec(true_vars, rf_cv_vars)
# 	rep_results[16, , nr] = c(obj$precision, obj$recall) 	
# 	obj = calc_prec_rec(true_vars, rf_point_vars)
# 	rep_results[17, , nr] = c(obj$precision, obj$recall)
# 	obj = calc_prec_rec(true_vars, rf_simul_vars)
# 	rep_results[18, , nr] = c(obj$precision, obj$recall)
#   
#   ##dynatree
# 	obj = calc_prec_rec(true_vars, dynatree_vars)
# 	rep_results[19, , nr] = c(obj$precision, obj$recall)
# 	
# 	##spikeslab 
# 	obj = calc_prec_rec(true_vars, spikeslab_vars)
# 	rep_results[20, , nr] = c(obj$precision, obj$recall)  
  
# 	write.csv(rep_results[, , nr], file = paste("partial_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, "_nr_", nr, "_revised.csv", sep = ""))	
	write.csv(rep_results[, , nr], file = paste("partial_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, "_nr_", nr, "_bart_gmax.csv", sep = ""))
}

results = matrix(0, nrow = NUM_ALGOS, ncol = 2)
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
		"RF_simul",
    "dynaTree",
    "spikeSlab") 
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
# write.csv(results, file = paste("complete_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, "_revised.csv", sep = ""))

###load results and print them to xtable
#sink("all_results_linear.tex")

#tryCatch(library(xtable), error = function(e){install.packages("xtable")}, finally = library(xtable))
#for (p in ps){
#	for (po_prop in po_props){
#		for (sigsq in sigsqs){			
#			p0 = ceiling(p * po_prop)
#			results = read.csv(paste("var_sel_sim_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
#			print(xtable(results, caption = paste("$p = ", p, ",~p_0 = ", p0, ",~\\sigsq = ", sigsq, "$\n", sep ="")))
#		}			
#	}
#}
#sink()
