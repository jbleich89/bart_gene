LAST_NAME = "kapelner"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/workspace/bart_gene/analysis_with_tuning/", sep = ""))
} else {
	setwd("analysis_with_tuning")
}

source("simulation_params.R")

#read in arguments supplied by qsub - this will tell use which gene to process
args = commandArgs(TRUE)
print(paste("args:", args))

if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}


#now load up the BART stuff
if (NOT_ON_GRID){
	setwd("C:/Users/Kapelner/workspace/CGMBART_GPL/")
} else {
	setwd("../CGMBART_GPL/")
}

source("r_scripts/bart_package.R")



find_important_tfs_for_gene_via_null_bart_sampling = function(gene, c_param, m_param){
	
	#pull the (adjusted) probabilities that the biologists give us and use as a prior
	cov_probs = prior_mat_adj[gene,]
	#as we discussed, we're going to weight the biologists guesses and bootstrap over them
	cov_prior_vec = 1 + as.numeric(c_param) * cov_probs #cov_prior = cov_prior / sum(cov_prior)	
	#pull out the responses and add them to the training data
	y = gene_train[, gene]
	training_data = data.frame(tf_train, y)
	#run BART many times and average the importances
	cat("get_averaged_true_var_props\n")
	var_props_avg = get_averaged_true_var_props(training_data, cov_prior_vec, m_param)
#	barplot(var_props_avg, names = names(var_props_avg), las = 2, main = paste("important tfs for gene", gene, "m =", m, "c =", c_param))
	
	#do the null permute sampling over many cores
	permute_mat = matrix(NA, nrow = NUM_PERMUTE_SAMPLES, ncol = ncol(tf_train))
	colnames(permute_mat) = colnames(tf_train)
	
	cat("get_null_permute_var_importances\n")
	for (b in 1 : NUM_PERMUTE_SAMPLES){
		permute_mat[b, ] = get_null_permute_var_importances(training_data, m_param)
	}
	
	
	#find tfs likely to be important - there are three methods to do this and a bunch of alpha cutoffs to choose from
	gene_result_c_m_methods_alphas = list()
	for (alpha in alphas){
		cat(paste("c:", c_param, "m:", m_param, "alpha:", alpha, "finding important tfs for gene:", gene, "\n"))
		
		gene_result_c_m_methods_alphas[[alpha]] = list()
		
		alpha_val = as.numeric(alpha)
		#pointwise
		pointwise_cutoffs = apply(permute_mat, 2, quantile, probs = 1 - alpha_val)
		important_tfs_at_alpha_pointwise = var_props_avg[var_props_avg > pointwise_cutoffs]
		
		#simult. max
		maxcut = quantile(apply(permute_mat, 1 ,max), 1 - alpha_val)
		important_tfs_at_alpha_simul_max = var_props_avg[var_props_avg >= maxcut]
		
		#simult. se
		perm_se = apply(permute_mat ,2 , sd)
		perm_mean = apply(permute_mat ,2 , mean)
		cover_constant = bisectK(tol = .01 , coverage = 1 - alpha_val, permute_mat = permute_mat, x_left = 1, x_right = 20, countLimit = 100, perm_mean = perm_mean, perm_se = perm_se)
		important_tfs_at_alpha_simul_se = var_props_avg[which(var_props_avg >= perm_mean + cover_constant * perm_se)]
		
		#now save them correctly
		gene_result_c_m_methods_alphas[[alpha]][["important_tfs_at_alpha_pointwise"]] = important_tfs_at_alpha_pointwise
		gene_result_c_m_methods_alphas[[alpha]][["important_tfs_at_alpha_simul_max"]] = important_tfs_at_alpha_simul_max
		gene_result_c_m_methods_alphas[[alpha]][["important_tfs_at_alpha_simul_se"]] = important_tfs_at_alpha_simul_se
	}
	
	gene_result_c_m_methods_alphas
}

get_averaged_true_var_props = function(training_data, cov_prior_vec, m){
	var_props = rep(0, ncol(training_data) - 1)
	for (i in 1 : NUM_REP_FOR_TRAIN){
		bart_machine = build_bart_machine(training_data, 
			num_trees = as.numeric(m), 
			num_burn_in = NUM_BURN_IN, 
			num_iterations_after_burn_in = NUM_ITER_AFTER, 
			cov_prior_vec = cov_prior_vec,
			run_in_sample = FALSE,
			verbose = FALSE)
		var_props = var_props + get_var_props_over_chain(bart_machine)
		destroy_bart_machine(bart_machine)
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
		num_trees = as.numeric(m), 
		num_burn_in = NUM_BURN_IN, 
		num_iterations_after_burn_in = NUM_ITER_AFTER,
		run_in_sample = FALSE,
		verbose = FALSE)
#	#just return the variable proportions	
	var_props = get_var_props_over_chain(bart_machine)
	destroy_bart_machine(bart_machine)
	var_props
}


#run for all genes and all c's and all m's
run_simulation_for_gene_cs_and_ms = function(gene_num){
	results = list()
	gene = gene_names[gene_num]
	results[[gene]] = list()
	for (c_param in cs){
		results[[gene]][[c_param]] = list()
		for (m_param in ms){
			results[[gene]][[c_param]][[m_param]] = list()
			
			cat(paste("c:", c_param, "m:", m_param, "finding important tfs for gene:", gene, "\n"))
			
			#now actually find the important tf's
			t0 = Sys.time()
			results[[gene]][[c_param]][[m_param]] = find_important_tfs_for_gene_via_null_bart_sampling(gene, c = c_param, m = m_param)
			
			cat(paste(" (completed in", round(Sys.time() - t0, 2), "min)\n"))
			save(results, file = paste("gene_results_", gene_num, ".RData", sep = ""))
		}
	}
	results
}



#go through results and validate

validate_for_gene_cs_and_ms = function(gene_num){
	gene = gene_names[gene_num]
	validation_oos_rmses = list()
	validation_oos_rmses[[gene]] = list()
	
	t0 = Sys.time()
	
	for (c_param in cs){		
		validation_oos_rmses[[gene]][[c_param]] = list()
		
		for (m_param in ms){
			validation_oos_rmses[[gene]][[c_param]][[m_param]] = list()
			
			for (method in METHODS){
				for (alpha in alphas){
					cat("validate for gene:", gene, "c:", c_param, "m:", m_param, "method:", method, "alpha:", alpha, "\n")
					#get important tfs for this combination of c, m, gene, and method
					important_tfs = names(results[[gene]][[c_param]][[m_param]][[alpha]][[method]])
					
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
						destroy_bart_machine(bart_machine)
						validation_oos_rmses[[gene]][[c_param]][[m_param]][[alpha]][[method]] = predict_obj$rmse						
					} else {
						L2_err = sum((y_cv - mean(y_train))^2)
						rmse = sqrt(L2_err / length(y_cv))
						validation_oos_rmses[[gene]][[c_param]][[m_param]][[alpha]][[method]] = rmse
					}					
				}				
			}
		}
		save(validation_oos_rmses, file = paste("validation_results_", GENE_NUM, ".RData", sep = ""))
	}
	cat(paste(" (completed in", round(Sys.time() - t0, 2), "min)\n"))
	
	validation_oos_rmses
}

if (NOT_ON_GRID){
	GENE_NUM = 770
}

#now run the simulation and validate
results = run_simulation_for_gene_cs_and_ms(GENE_NUM)
#load(paste("gene_results_", GENE_NUM, ".RData", sep = ""))
validation_oos_rmses = validate_for_gene_cs_and_ms(GENE_NUM)



