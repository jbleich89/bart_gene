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
ps = c(10, 20, 100, 200, 1000)
po_props = c(0.1, 0.2, 0.5, 0.75, 1)
sigsqs = c(0.1, 0.5, 1, 5, 10)

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
	set_bart_num_cores(4)
}

p = param_mat[iter_num, 1]
po_prop = param_mat[iter_num, 2]
sigsq = param_mat[iter_num, 3]


#generate data
X = matrix(runif(n * p), ncol = p)
p0 = round(p * po_prop)
beta_vec = c(rep(1, p0), rep(0, p - p0))
error = rnorm(n, 0, sqrt(sigsq))
y = as.numeric(X %*% beta_vec + error)


#now build bart machine (even though we don't really have to, but it's nice to have)

bart_machine = build_bart_machine(as.data.frame(X), y)






#y = 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5] + 
#		10 * sin(pi * X[, 6] * X[, 7]) + 20 * (X[, 8] - 0.5)^2 + 10 * X[, 9] + 5 * X[, 10] +
#		error