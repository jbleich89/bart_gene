if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/workspace/bart_gene/analysis_with_tuning", sep = ""))
} else {
	setwd("../bart_gene/analysis_with_tuning")
}

source("simulation_params.R")

n_test = nrow(tf_test)
#make sse matrix

all_sse_results = n_test * all_rmse_results^2
sse_reductions_per_var = matrix(NA, nrow = nrow(all_rmse_results), ncol = ncol(all_rmse_results))
rownames(sse_reductions_per_var) = rownames(all_rmse_results)
colnames(sse_reductions_per_var) = colnames(all_rmse_results)

for (g in 1 : nrow(all_sse_results)){
	gene_name = rownames(all_sse_results)[g]
	gene_test_data = gene_test[, gene_name]
	sse_null = sum((gene_test_data - mean(gene_test_data))^2)
	sse_reductions_per_var[g, ] = (sse_null - all_sse_results[g, ]) / all_num_var_results_no_na[g, ]
}

head(sse_reductions_per_var)
head(all_rmse_results)

boxplot(sse_reductions_per_var, ylim = c(-10, 10), ylab = "out-of-sample Scaled SSE reduction", main = "Out-of-Sample Scaled SSE reduction by Method")
points(apply(all_rmse_results, 2, mean), pch = "-", col = "blue", cex = 5)
abline(a = mean(all_rmse_results[, 1]), b = 0, col = "red")





##### rmse per num vars

all_rmse_per_num_var = all_rmse_results / all_num_var_results_no_na
#
boxplot(all_rmse_per_num_var, ylab = "out-of-sample Scaled SSE reduction", main = "Out-of-Sample Scaled SSE reduction by Method", ylim = c(0, 0.5))
points(apply(all_rmse_per_num_var, 2, mean), pch = "-", col = "blue", cex = 5)
abline(a = mean(all_rmse_per_num_var[, 1]), b = 0, col = "red")



all_rmse_minus_null_per_num_var = (all_rmse_results - matrix(rep(all_rmse_results[, 1], ncol(all_rmse_results)), ncol = ncol(all_rmse_results))) / all_num_var_results_no_na
#
boxplot(all_rmse_minus_null_per_num_var, ylab = "out-of-sample Scaled SSE reduction", main = "Out-of-Sample Scaled SSE reduction by Method")
points(apply(all_rmse_minus_null_per_num_var, 2, mean), pch = "-", col = "blue", cex = 5)
abline(a = mean(all_rmse_minus_null_per_num_var[, 1]), b = 0, col = "red")