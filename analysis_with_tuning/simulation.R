if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/workspace/bart_gene/")
}

source("analysis_with_tuning/simulation_params.R")

#read in arguments supplied by qsub - this will tell use which gene to process
args = commandArgs(TRUE)
print(paste("args:", args))

if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}

#get necessary functions and data


source("helper_functions.R")

if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene")
} else {
	setwd("../shane_data")
}

priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

priorMat = as.matrix(priors[, 3 : ncol(priors)])
rownames(priorMat) = priors[, 1] 
colnames(priorMat) = colnames(priors)[3 : ncol(priors)]
prior_mat_adj = prior_adj(priorMat)
result = setup() ##5 objects
gene_train = result[["gene.train"]]
gene_cv = result[["gene.cv"]]
gene_test = result[["gene.test"]]##gene.train and gene.test rows are obs and cols are genes
tf_train = result[["tf.train"]]
tf_cv = result[["tf.cv"]]
tf_test = result[["tf.test"]]##TF train and TF test: rows are obs and cols are TFs
gene_names = as.character(gene.exp[, 2]) ##gene names

#now load up the BART stuff
if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/workspace/CGMBART_GPL/")
} else {
	setwd("../CGMBART_GPL/")
}

source("r_scripts/bart_package.R")



find_important_tfs_for_gene_via_null_bart_sampling = function(gene, c_param, m, ...){
	t0 = Sys.time()

	#pull the (adjusted) probabilities that the biologists give us and use as a prior
	cov_probs = prior_mat_adj[gene,]
	#as we discussed, we're going to weight the biologists guesses and bootstrap over them
	cov_prior_vec = 1 + c_param * cov_probs #cov_prior = cov_prior / sum(cov_prior)	
	#pull out the responses and add them to the training data
	y = gene_train[, gene]
	training_data = data.frame(tf_train, y)
	#run BART many times and average the importances
	var_props_avg = get_averaged_true_var_props(training_data, cov_prior_vec, m)
	barplot(var_props_avg, names = names(var_props_avg), las = 2, main = paste("important tfs for gene", gene, "m =", m, "c =", c_param))
	
	#do the null permute sampling over many cores
	permute_mat = matrix(NA, nrow = NUM_PERMUTE_SAMPLES, ncol = ncol(tf_train))
	colnames(permute_mat) = colnames(tf_train)
	
	for (b in 1 : NUM_PERMUTE_SAMPLES){
		permute_mat[b, ] = get_null_permute_var_importances(training_data, m)
	}
	
	
	#find tfs likely to be important - there are three methods to do this
	
	#pointwise
	pointwise_cutoffs = apply(permute_mat, 2, quantile, probs = 1 - ALPHA)
	important_tfs_at_alpha_pointwise = var_props_avg[var_props_avg > pointwise_cutoffs]
	
	#simult. max
	maxcut = quantile(apply(permute_mat, 1 ,max), 1 - ALPHA)
	important_tfs_at_alpha_simul_max = var_props_avg[var_props_avg >= maxcut]
	
	#simult. se
	perm_se = apply(permute_mat ,2 , sd)
	perm_mean = apply(permute_mat ,2 , mean)
	cover_constant = bisectK(tol = .01 , coverage = 1- ALPHA, permute_mat = permute_mat, x_left = 1, x_right = 20, countLimit = 100, perm_mean = perm_mean, perm_se = perm_se)
	important_tfs_at_alpha_simul_se = var_props_avg[which(var_props_avg >= perm_mean + cover_constant * perm_se)]
	
	list(
		gene = gene,
		var_props_avg = round(var_props_avg, 3),
		important_tfs_at_alpha_pointwise = important_tfs_at_alpha_pointwise,
		important_tfs_at_alpha_simul_max = important_tfs_at_alpha_simul_max,
		important_tfs_at_alpha_simul_se = important_tfs_at_alpha_simul_se,
		time_elapsed = round(Sys.time() - t0, 2)
	)
}

get_averaged_true_var_props = function(training_data, cov_prior_vec, m){
	var_props = rep(0, ncol(training_data) - 1)
	for (i in 1 : NUM_REP_FOR_TRAIN){
		bart_machine = build_bart_machine(training_data, 
			num_trees = m, 
			num_burn_in = NUM_BURN_IN, 
			num_iterations_after_burn_in = NUM_ITER_AFTER, 
			cov_prior_vec = cov_prior_vec,
			run_in_sample = FALSE,
			verbose = FALSE)
		var_props = var_props + get_var_props_over_chain(bart_machine)
#		plot_sigsqs_convergence_diagnostics(bart_machine)
	}
	#average over many runs
	var_props / NUM_REP_FOR_TRAIN
}

get_null_permute_var_importances = function(training_data, m){
	#permute the responses to disconnect x and y
	training_data$y = sample(training_data$y, replace = FALSE)
	#build BART on this permuted training data
	bart_machine = build_bart_machine(training_data, 
		num_trees = m, 
		num_burn_in = NUM_BURN_IN, 
		num_iterations_after_burn_in = NUM_ITER_AFTER,
		run_in_sample = FALSE,
		verbose = FALSE)
#	#just return the variable proportions	
	get_var_props_over_chain(bart_machine)
}


#run for all genes and all c's and all m's
run_simulation_for_all_genes_cs_and_ms = function(gene_num){
	results = list()
	gene = gene_names[gene_num]
	results[[gene]] = list()
	for (c_param in cs){ 
		results[[gene]][[as.character(c_param)]] = list()
		for (m_param in ms){
			results[[gene]][[as.character(c_param)]][[as.character(m_param)]] = list()
			cat(paste("c:", c_param, "m:", m_param, "finding important tfs for gene:", gene, "\n"))
			#now actually find the important tf's
			result = find_important_tfs_for_gene_via_null_bart_sampling(gene, c = c_param, m = m_param)
			
			for (method in METHODS){
				results[[gene]][[as.character(c_param)]][[as.character(m_param)]][[method]] = result[[method]]
			}
			cat(paste(" (completed in", result[["time_elapsed"]], "min)\n"))
			save(results, file = paste("gene_results_", gene_num, ".RData", sep = ""))
		}
	}
	results
}



#go through results and validate

validate_for_all_genes_cs_and_ms = function(gene_num){
	gene = gene_names[gene_num]
	validation_oos_rmses = list()
	validation_oos_rmses[[gene]] = list()
	
	for (c_param in cs){		
		validation_oos_rmses[[gene]][[as.character(c_param)]] = list()
		
		for (m_param in ms){
			validation_oos_rmses[[gene]][[as.character(c_param)]][[as.character(m_param)]] = list()
			
			for (method in METHODS){
				#get important tfs for this combination of c, m, gene, and method
				important_tfs = names(results[[gene]][[as.character(c_param)]][[as.character(m_param)]][[method]])
				
				y_train = gene_train[, gene]
				y_cv = gene_cv[, gene]
				
				if (length(important_tfs) > 0){
					#build training data
					
					training_data = data.frame(tf_train, y = y_train)
					training_data = training_data[, c(important_tfs, "y")]
					
					#now run BART model with 200 trees
					bart_machine = build_bart_machine(training_data, 
						num_trees = NUM_TREES_FOR_EVAL, 
						num_burn_in = NUM_BURN_IN, 
						num_iterations_after_burn_in = NUM_ITER_AFTER,
						verbose = FALSE)
					
					#predict on cv data set only with important tf's						
					cv_data = data.frame(tf_cv, y = y_cv)
					cv_data = cv_data[, c(important_tfs, "y")]
					predict_obj = bart_predict_for_test_data(bart_machine, cv_data)
					validation_oos_rmses[[gene]][[as.character(c_param)]][[as.character(m_param)]][[method]] = predict_obj$rmse						
				} else {
					L2_err = sum((y_cv - mean(y_train))^2)
					rmse = sqrt(L2_err / length(y_train))						
					validation_oos_rmses[[gene]][[as.character(c_param)]][[as.character(m_param)]][[method]] = rmse
				}
			}
		}
		save(validation_oos_rmses, file = paste("validation_results_", GENE_NUM, ".RData", sep = ""))
	}
	validation_oos_rmses
}

#now run the simulation and validate
results = run_simulation_for_all_genes_cs_and_ms(GENE_NUM)
#load(paste("gene_results_", GENE_NUM, ".RData", sep = ""))
validation_oos_rmses = validate_for_all_genes_cs_and_ms(GENE_NUM)
