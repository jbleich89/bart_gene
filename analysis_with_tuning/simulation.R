tryCatch(library(parallel), error = function(e){install.packages("parallel")}, finally = library(parallel))

#simulation params

cs = c(
	0,
	0.5,
	1,
	1.5,
	2,
	3
)

m = c(
	10,
	20
)

NUM_PERMUTE_SAMPLES = 50
NUM_BURN_IN = 5000
NUM_ITER_AFTER = 2500
NUM_REP_FOR_TRAIN = 10
NUM_CORES = 1

#get necessary functions and data

setwd("C:/Users/Kapelner/workspace/bart_gene")
#source("real_functions.R")
source("helper_functions.R")

setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene")

priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

priorMat = as.matrix(priors[, 3 : ncol(priors)])
rownames(priorMat) = priors[,1] 
colnames(priorMat) = colnames(priors)[3 : ncol(priors)]
prior_mat_adj = prior_adj(priorMat)
out = setup() ##5 objects
gene_train = out[["gene.train"]]
gene_cv = out[["gene.cv"]]
gene_test = out[["gene.test"]]##gene.train and gene.test rows are obs and cols are genes
tf_train = out[["tf.train"]]
tf_cv = out[["tf.cv"]]
tf_test = out[["tf.test"]]##TF train and TF test: rows are obs and cols are TFs
gene_names = as.character(gene.exp[, 2]) ##gene names

#now load up the BART stuff
setwd("C:/Users/Kapelner/workspace/CGMBART_GPL/")
source("r_scripts/bart_package.R")

results = list()
for (gene in gene_names[1 : 100]){
	cat(paste("finding important tfs for gene:", gene))
	results[[gene]] = find_important_tfs_for_gene_via_null_bart_sampling(gene)
}


find_important_tfs_for_gene_via_null_bart_sampling = function(gene, ...){
	t0 = Sys.time()
	output = list()
	out[[gene]] = list()
	out[[gene]][["name"]] = gene
	cov_probs = prior_mat_adj[gene,]
	cov_prior = 1 + c * cov_probs 
	y = gene_train[, gene]
	training_data = data.frame(tf_train, y)
	var_props = get_averaged_true_var_props(training_data, cov_prior)
	barplot(var_props, names = names(var_props), las = 2, main = paste("important tfs for gene", gene))
	
	
	out[["var_props_avg"]] = round(var_props, 3)
	
#	if (runBoot){
		permute_list = mclapply(1 : NUM_PERMUTE_SAMPLES, get_permute_iter, training_data = training_data, mc.cores = NUM_CORES)
		permute_mat = t(do.call(cbind, permute_list))
		colnames(permute_mat) = colnames(tf_train)
		
		#Pointwise
		pointwise_cutoffs = apply(permute_list, 2, quantile, probs = .95)
		pointwise_vars = var_props[which(var_props > pointwise_cutoffs)]
		pointwise_vars
		##Simult. Max 
#	}
}

get_averaged_true_var_props = function(training_data, cov_prior){
	var_props = rep(0, ncol(training_data) - 1)
	for (i in 1 : NUM_REP_FOR_TRAIN){
		bart_machine = build_bart_machine(training_data, 
				num_trees = m, 
				num_burn_in = NUM_BURN_IN, 
				num_iterations_after_burn_in = NUM_ITER_AFTER, 
				cov_prior_vec = cov_prior,
				verbose = FALSE)
		var_props = var_props + get_var_props_over_chain(bart_machine)
#		plot_sigsqs_convergence_diagnostics(bart_machine)
	}
	#average over many runs
	var_props / Nsim
}

get_permute_iter = function(iter, training_data){
	#permute the responses to disconnect x and y
	training_data$y = sample(training_data$y, replace = FALSE)
	#build BART on this permuted training data
	bart_machine = build_bart_machine(training_data, 
			num_trees = m, 
			num_burn_in = NUM_BURN_IN, 
			num_iterations_after_burn_in = NUM_ITER_AFTER,
			verbose = FALSE)
	#just return the variable proportions	
	get_var_props_over_chain(bart_machine)
}














#gene = gene_names[[3]]
#y = gene.train[, gene]
#training_data = data.frame(tf.train, y)
#
#graphics.off()
#for (j in 1 : 3){
#	Nsim = 10
#	var_props = rep(0, ncol(training_data) - 1)
#	for (i in 1 : Nsim){
#		bart_machine = build_bart_machine(training_data, 
#				num_trees = m, 
#				num_burn_in = NUM_BURN_IN, 
#				num_iterations_after_burn_in = NUM_ITER_AFTER, 
#				cov_prior_vec = cov_prior)
#		var_props = var_props + get_var_props_over_chain(bart_machine)
#	}
#	#average over many runs
#	var_props = var_props / Nsim
#	windows()
#	barplot(var_props, names = names(var_props), las = 2)
#}