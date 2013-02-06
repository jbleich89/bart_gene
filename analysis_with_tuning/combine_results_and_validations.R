#to be run on GRID only!!!

options(width = 150)
MAX_GENE_NUM = 1150

if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/workspace/bart_gene/")
}

source("analysis_with_tuning/simulation_params.R")


if (.Platform$OS.type == "windows"){
	setwd("C:\\Users\\kapelner\\workspace\\CGMBART_GPL")
} else {
	setwd("../CGMBART_GPL")
}


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

save(all_results, file = "all_results.RData")
save(all_validations, file = "all_validations.RData")

