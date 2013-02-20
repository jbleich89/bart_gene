####win by THREHSOLD% matrix

THRESHOLD = 0.2

wins = list()

for (gene_num in 1 : length(all_results)){
	gene_name = names(all_results)[gene_num]
	wins[[gene_name]] = array(0, c(length(cs), length(METHODS), length(alphas)))
	dimnames(wins[[gene_name]]) = list(
			"c" = as.character(cs),
			"method" = METHODS,
			"alpha" = alphas
	) 
	
	best = 9999999
	second_best = 9999999
	lowest_c = NULL
	lowest_method = NULL
	lowest_alpha = NULL
	
	all_rmses = as.data.frame(matrix(NA, nrow = length(c_param) * length(METHODS) * length(alphas), ncol = 3))
	names(all_rmses) = c("rmse", "c", "m")
	counter = 1
	for (c_param in cs){	
		c_param = as.character(c_param)
		for (method in METHODS){
			for (alpha in alphas){
				oo_rmse = all_validations[[gene_name]][[c_param]][["20"]][[alpha]][[method]]
				all_rmses[counter, ] = c(oo_rmse, c_param, method)
#				if (length(oo_rmse) == 0){
#					cat(paste("ERROR:", gene_name, c_param, method, "\n"))
#				} else if (oo_rmse < best){
#					
#					second_best = best
#					best = oo_rmse
#					lowest_c = c_param
#					lowest_method = method
#					lowest_alpha = alpha
#					print(paste(gene_num, "oo_rmse beaten  best:", best, "second_best", second_best, "c", c_param, "method", method))
#				}
				counter = counter + 1
			}
		}
	}	
	
	#sort
	all_rmses$rmse = as.numeric(all_rmses$rmse) 
	all_rmses = all_rmses[sort.int(all_rmses$rmse, index.return = TRUE)$ix, ]
	
	if ((all_rmses$rmse[2] - all_rmses$rmse[1]) / all_rmses$rmse[2] >= THRESHOLD){
#		print(paste(gene_num, "FINAL oo_rmse MORE than THRESHOLD best:", all_rmses$rmse[1], "second_best", all_rmses$rmse[2], "c", all_rmses$c[1], "method", all_rmses$m[1]))
		wins[[gene_name]][all_rmses$c[1], all_rmses$m[1], alphas[1]] = 1
	}	
}

aggregated_win_by_threshold = array(0, c(length(cs),length(METHODS), length(alphas)))
dimnames(aggregated_win_by_threshold) = list(
		"c" = as.character(cs),
		"method" = METHODS,
		"alpha" = alphas
) 

for (gene_num in 1 : length(all_results)){
	gene_name = names(all_results)[gene_num]
	aggregated_win_by_threshold = aggregated_win_by_threshold + wins[[gene_name]]
}
aggregated_win_by_threshold
#num genes where it mattered
sum(aggregated_win_by_threshold)
#pct genes where it mattered
sum(aggregated_win_by_threshold) / length(all_results)

xtable(aggregated_win_by_threshold[,, 1], digits = 0)
