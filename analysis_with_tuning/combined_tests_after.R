

#run on the GRID ONLY!!
setwd("../CGMBART_GPL")

MAX_GENE_NUM = 4000

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
tail(all_rmse_results)
tail(all_num_var_results)

save(all_rmse_results, file = paste("all_rmse_results.RData", sep = ""))
save(all_num_var_results, file = paste("all_num_var_results.RData", sep = ""))


#########load results into R at this point ON THE LOCAL MACHINE

#lasso is worse
sum((all_rmse_results[, "Lasso"] - all_rmse_results[, "BART-Best"]) / all_rmse_results[, "BART-Best"] > 0.05)
#bart is worse
sum((all_rmse_results[, "BART-Best"] - all_rmse_results[, "Lasso"]) / all_rmse_results[, "BART-Best"] > 0.05)


sum(all_rmse_results[, "BART-Best"] < all_rmse_results[, "Lasso"])

#oos rmses for all methods
boxplot(all_rmse_results, ylim = c(0, 1.1), ylab = "out-of-sample RMSE", main = "Out-of-Sample RMSE by Method")
points(apply(all_rmse_results, 2, mean), pch = "-", col = "blue", cex = 5)
abline(a = mean(all_rmse_results[, 1]), b = 0, col = "red")

all_num_var_results_no_na = t(sapply(1 : nrow(all_num_var_results), function(i){ifelse(is.na(all_num_var_results[i, ]), 39, all_num_var_results[i, ])}))

boxplot(all_num_var_results_no_na, ylab = "Num TF's Selected", main = "Num TF's Selected by Method")
points(apply(all_num_var_results_no_na, 2, mean), pch = "-", col = "blue", cex = 5)
