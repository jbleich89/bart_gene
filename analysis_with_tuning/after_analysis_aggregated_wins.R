#now we want to see who won...
wins = list()

for (gene_num in 1 : length(all_results)){
	gene_name = names(all_results)[gene_num]
	wins[[gene_name]] = array(0, c(length(cs), length(METHODS), length(alphas)))
	dimnames(wins[[gene_name]]) = list(
			"c" = as.character(cs),
			"method" = METHODS,
			"alpha" = alphas
	) 
	
	min = 9999999
	lowest_c = NULL
	lowest_method = NULL
	lowest_alpha = NULL
	for (c_param in cs){	
		c_param = as.character(c_param)
		for (method in METHODS){
			for (alpha in alphas){
				oo_rmse = all_validations[[gene_name]][[c_param]][["20"]][[alpha]][[method]]
				if (length(oo_rmse) == 0){
					cat(paste("ERROR:", gene_name, c_param, method, "\n"))
				} else if (oo_rmse < min){
					min = oo_rmse
					lowest_c = c_param
					lowest_method = method
					lowest_alpha = alpha
				}
			}
		}
	}	
	wins[[gene_name]][lowest_c, lowest_method, lowest_alpha] = 1
}

aggregated_wins = array(0, c(length(cs),length(METHODS), length(alphas)))
dimnames(aggregated_wins) = list(
		"c" = as.character(cs),
		"method" = METHODS,
		"alpha" = alphas
) 

for (gene_num in 1 : length(all_results)){
	gene_name = names(all_results)[gene_num]
	aggregated_wins = aggregated_wins + wins[[gene_name]]
}
aggregated_wins
sum(aggregated_wins)

##sum aggregated wins

##prior stuff 
library(xtable)
x = aggregated_wins[,1,] + aggregated_wins[,2,]+ aggregated_wins[,3,]
xtable(prop.table(as.table(x)))
