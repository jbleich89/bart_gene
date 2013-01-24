

if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/workspace/bart_gene/")
}

source("analysis_with_tuning/simulation_params.R")


if (.Platform$OS.type == "windows"){
	setwd("C:\\Users\\kapelner\\workspace\\CGMBART_GPL")
} else {
	setwd("../CGMBART_GPL")
}


MAX_GENE_NUM = 100

all_results = list()
all_validations = list()

for (gene_num in 1 : MAX_GENE_NUM){
	#get the results first
	load(paste("gene_results_", gene_num, ".RData", sep = ""))
	gene_name = names(results)[1]
	all_results[[gene_name]] = results[[gene_name]]
	
	#get the validations second
	load(paste("validation_results_", gene_num, ".RData", sep = ""))
	all_validations[[gene_name]] = validation_oos_rmses[[gene_name]]
}

#now we want to see who won...

wins = list()

for (gene_num in 1 : length(all_results)){
	gene_name = names(all_results)[gene_num]
	wins[[gene_name]] = matrix(0, nrow = length(cs), ncol = length(METHODS))
	rownames(wins[[gene_name]]) = as.character(cs)
	colnames(wins[[gene_name]]) = METHODS
	
	min = 9999999
	lowest_c = NULL
	lowest_method = NULL
	for (c_param in c(0,1,2)){	
		c_param = as.character(c_param)
		for (method in METHODS){
			oo_rmse = all_validations[[gene_name]][[c_param]][["20"]][[method]]
			if (oo_rmse < min){
				min = oo_rmse
				lowest_c = c_param
				lowest_method = method
			}
		}
	}	
	wins[[gene_name]][lowest_c, lowest_method] = 1
}

aggregated_wins = matrix(0, nrow = length(cs), ncol = length(METHODS))
rownames(aggregated_wins) = as.character(cs)
colnames(aggregated_wins) = METHODS

for (gene_num in 1 : length(all_results)){
	gene_name = names(all_results)[gene_num]
	aggregated_wins = aggregated_wins + wins[[gene_name]]
}