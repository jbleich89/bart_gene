

#run on the GRID ONLY!!
setwd("../CGMBART_GPL")

MAX_GENE_NUM = 1000

simulation_names = c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso-CV", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full", "Rob-Best")
all_rmse_results = matrix(NA, nrow = MAX_GENE_NUM, ncol = length(simulation_names))
colnames(all_rmse_results) = simulation_names
rownames(all_rmse_results) = 1 : MAX_GENE_NUM
all_num_var_results = matrix(NA, nrow = MAX_GENE_NUM, ncol = length(simulation_names))
colnames(all_num_var_results) = simulation_names
rownames(all_num_var_results) = 1 : MAX_GENE_NUM


for (g in 1 : MAX_GENE_NUM){
	load(file = paste("rmse_results_", g, ".RData", sep = ""))
	load(file = paste("num_var_results_", g, ".RData", sep = ""))
	all_rmse_results[g, ] = rmse_results
	rownames(all_rmse_results)[g] = rownames(rmse_results)
	all_num_var_results[g, ] = num_var_results
	rownames(all_num_var_results)[g] = rownames(num_var_results)[1]
}

head(all_rmse_results)
head(all_num_var_results)

save(all_rmse_results, file = paste("all_rmse_results.RData", sep = ""))
save(all_num_var_results, file = paste("all_num_var_results.RData", sep = ""))

sum((all_rmse_results[, "Lasso"] - all_rmse_results[, "BART-Best"]) / all_rmse_results[, "BART-Best"] > 0.05)

sum(all_rmse_results[, "BART-Best"] < all_rmse_results[, "Lasso"])
