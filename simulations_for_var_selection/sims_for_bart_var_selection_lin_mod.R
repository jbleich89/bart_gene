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

#if (NOT_ON_GRID){
#	setwd(paste("C:/Users/", LAST_NAME, "/workspace/bart_gene/simulations_for_var_selection/", sep = ""))
#} else {
#	setwd("simulations_for_var_selection")
#}

###
n = 250
ps = c(10, 20, 100, 200)
po_props = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.75)
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


#generate linear model data
X = matrix(runif(n * p), ncol = p)
p0 = ceiling(p * po_prop)
true_vars = 1 : p0
beta_vec = c(rep(1, p0), rep(0, p - p0))
error = rnorm(n, 0, sqrt(sigsq))
y = as.numeric(X %*% beta_vec + error)
X = as.data.frame(X)
colnames(X) = seq(1, p)

sqrt(sigsq) / (sum(abs(as.matrix(X) %*% beta_vec)) / n)


#now build bart machine (even though we don't really have to, but it's nice to have)

bart_machine = build_bart_machine(X, y, num_trees = 1)

#do var selection with bart
bart_variables_select_obj = var_selection_by_permute_response(bart_machine, num_permute_samples = 100)
bart_ptwise_vars = sort(as.numeric(bart_variables_select_obj$important_vars_pointwise))
bart_simul_max_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_max))
bart_simul_se_vars = sort(as.numeric(bart_variables_select_obj$important_vars_simul_se))

#do var selection with stepwise
step_reg = stepAIC(object = lm(y ~ ., data = X), direction = "backward" , trace = F)
stepwise_vars = names(step_reg$coefficients)
stepwise_vars = as.numeric(gsub("`", "", stepwise_vars))
stepwise_vars = sort(stepwise_vars[!is.na(stepwise_vars)])

#do var selection with lasso
lasso_matrix_vars = coef(cv.glmnet(as.matrix(X), y, alpha = 1 , nfolds = 5, type.measure = "mse"))
lasso_matrix_vars = which(lasso_matrix_vars != 0) - 1
lasso_matrix_vars = sort(lasso_matrix_vars[lasso_matrix_vars != 0]) #kill intercept if it exists


#### now what did we get right?
results = matrix(NA, nrow = 5, ncol = 3)
rownames(results) = c("BART_pointwise", "BART_simul_max", "BART_simul_se", "stepwise", "lasso")
colnames(results) = c("precision", "recall", "F1_measure")

calc_prec_rec_F1 = function(true_vars, regression_vars){
	true_vars_found = intersect(true_vars, regression_vars)
	tps = length(true_vars_found)
	fps = length(setdiff(regression_vars, true_vars_found))
	fns = length(setdiff(true_vars, true_vars_found))
	precision = tps / (tps + fps)
	recall = tps / (tps + fns)
	list(
		precision = precision,
		recall = recall,
		F1_measure = 2 * precision * recall / (precision + recall)
	)
}

obj = calc_prec_rec_F1(true_vars, bart_ptwise_vars)
results[1, ] = c(obj$precision, obj$recall, obj$F1_measure)
obj = calc_prec_rec_F1(true_vars, bart_simul_max_vars)
results[2, ] = c(obj$precision, obj$recall, obj$F1_measure)
obj = calc_prec_rec_F1(true_vars, bart_simul_se_vars)
results[3, ] = c(obj$precision, obj$recall, obj$F1_measure)
obj = calc_prec_rec_F1(true_vars, stepwise_vars)
results[4, ] = c(obj$precision, obj$recall, obj$F1_measure)
obj = calc_prec_rec_F1(true_vars, lasso_matrix_vars)
results[5, ] = c(obj$precision, obj$recall, obj$F1_measure)

#save results
write.csv(results, file = paste("../bart_gene/simulations_for_var_selection/var_sel_sim_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))

###load results and print them to xtable
#sink("all_results.tex")
#tryCatch(library(xtable), error = function(e){install.packages("xtable")}, finally = library(xtable))
#for (p in ps){
#	for (po_prop in po_props){
#		for (sigsq in sigsqs){
#			if (p != 1000){				
#				p0 = round(p * po_prop)
#				results = read.csv(paste("var_sel_sim_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
#				print(xtable(results, caption = paste("$p = ", p, ", p_0 = ", p0, ", \\sigsq = ", sigsq, "$\n", sep ="")))
#			}
#		}			
#	}
#}
#sink()

#y = 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5] + 
#		10 * sin(pi * X[, 6] * X[, 7]) + 20 * (X[, 8] - 0.5)^2 + 10 * X[, 9] + 5 * X[, 10] +
#		error