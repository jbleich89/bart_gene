##BART tests


##Needs to take gene number as an argument
gene_num = 1

setwd("~/Research_Genomics/")

priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

source("~/Research_Genomics/bart_gene/analysis_with_tuning/simulation_params_JB.R")


setwd("C://Users/jbleich/Desktop/Dropbox/BART_gene/")

load("all_results.RData")
load("all_validations.RData")

setwd("C://Users/jbleich/workspace/CGMBART_GPL/")
source("r_scripts/bart_package.R")


rmse_list = all_validations[[gene_num]]
validation_list = all_results[[gene_num]]
rm(all_results) ; rm(all_validations)


min_rmse = Inf
tf_list = c()
for(c in cs){
  for(method in METHODS){
    if(rmse_list[[as.character(c)]][["20"]][[method]] < min_rmse){
      best_tfs = names(validation_list[[as.character(c)]][["20"]][[method]])
      min_rmse = rmse_list[[as.character(c)]][["20"]][[method]]
    }
  }
}

y_train = gene_train[, gene_num]
y_test = gene_test[, gene_num]

if (length(best_tfs) > 0){
  #build training data

  training_data = data.frame(tf_train, y = y_train)
  training_data = training_data[, c(best_tfs, "y")]
  
  #now run BART model with 200 trees
  bart_machine = build_bart_machine(training_data, 
                                    num_trees = NUM_TREES_FOR_EVAL, 
                                    num_burn_in = NUM_BURN_IN, 
                                    num_iterations_after_burn_in = NUM_ITER_AFTER, 
                                    verbose = FALSE)
  
  #predict on cv data set only with important tf's						
  test_data = data.frame(tf_test, y = y_test)
  test_data = test_data[, c(best_tfs, "y")]
  predict_obj = bart_predict_for_test_data(bart_machine, test_data)
  rmse = predict_obj$rmse						
  rmse
} else {
  L2_err = sum((y_test - mean(y_train))^2)
  rmse = sqrt(L2_err / length(y_train))						
}



