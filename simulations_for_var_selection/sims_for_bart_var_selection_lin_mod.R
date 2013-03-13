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
	list(
		precision = tps / (tps + fps),
		recall = tps / (tps + fns)
	)
}

###
num_replicates = 10
n = 250
ps = c(20, 100, 200, 500, 1000)
po_props = c(0.01, 0.05, 0.1, 0.2)
sigsqs = c(0.1, 0.5, 1, 5)

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
po_prop = param_mat[iter_num, 2]
sigsq = param_mat[iter_num, 3]


rep_results = array(NA, c(6, 2, num_replicates))

######replicate a few times
for (nr in 1 : num_replicates){
	#generate linear model data
	X = matrix(rnorm(n * p), ncol = p)
	p0 = ceiling(p * po_prop)
	true_vars = 1 : p0
	beta_vec = c(rep(1, p0), rep(0, p - p0))
	error = rnorm(n, 0, sqrt(sigsq))
	y = as.numeric(X %*% beta_vec + error)
	X = as.data.frame(X)
	colnames(X) = seq(1, p)
	
#	sqrt(sigsq) / (sum(abs(as.matrix(X) %*% beta_vec)) / n)	
	
	#now build bart machine (even though we don't really have to, but it's nice to have)
	
	bart_machine = build_bart_machine(X, y, num_trees = 1, run_in_sample = FALSE)
	
	#do var selection with bart
	bart_variables_select_obj = var_selection_by_permute_response(bart_machine, plot = FALSE)
	bart_ptwise_vars = sort(as.numeric(bart_variables_select_obj$important_vars_pointwise))
	bart_simul_max_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_max))
	bart_simul_se_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_se))
	
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
	obj = calc_prec_rec(true_vars, bart_ptwise_vars)
	rep_results[1, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_max_vars)
	rep_results[2, , nr] = c(obj$precision, obj$recall)
	obj = calc_prec_rec(true_vars, bart_simul_se_vars)
	rep_results[3, , nr] = c(obj$precision, obj$recall)
	if (p < n){
		obj = calc_prec_rec(true_vars, stepwise_backward_vars)
		rep_results[4, , nr] = c(obj$precision, obj$recall)
		obj = calc_prec_rec(true_vars, stepwise_forward_vars)
		rep_results[5, , nr] = c(obj$precision, obj$recall)		
	}
	obj = calc_prec_rec(true_vars, lasso_matrix_vars)
	rep_results[6, , nr ] = c(obj$precision, obj$recall)
}

results = matrix(0, nrow = 6, ncol = 2)
rownames(results) = c("BART_pointwise", "BART_simul_max", "BART_simul_se", "stepwise_backward",  "stepwise_forward", "lasso")
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
write.csv(results, file = paste("../bart_gene/simulations_for_var_selection/var_sel_sim_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))

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
