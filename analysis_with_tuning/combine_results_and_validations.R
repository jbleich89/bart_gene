options(width = 150)
MAX_GENE_NUM = 999

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

save(all_validations, file = "all_validations.RData")

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
	for (c_param in cs){	
		c_param = as.character(c_param)
		for (method in METHODS){
			
			oo_rmse = all_validations[[gene_name]][[c_param]][["20"]][[method]]
			if (length(oo_rmse) == 0){
				cat(paste("ERROR:", gene_name, c_param, method, "\n"))
			} else if (oo_rmse < min){
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
aggregated_wins
sum(aggregated_wins)


#now we want a boxplot

rmses = array(NA, c(length(cs), length(METHODS), MAX_GENE_NUM))

for (gene_num in 1 : MAX_GENE_NUM){
	validations = all_validations[gene_num]
	for (i in 1 : length(cs)){	
		c_param = as.character(cs[i])
		for (j in 1 : length(METHODS)){
			rmses[i, j, ] = sapply(1 : length(all_validations), function(s){all_validations[[s]][[c_param]][["20"]][[METHODS[j]]]})
		}
	}	
	
}

names = c()
for (c_param in cs){
	for (method in METHODS){
		names = c(names, paste(method, "c =", c_param))
	}
}

par(mar = c(20,3,3,3))
boxplot(rmses[1, 1, ], rmses[1, 2, ], rmses[1, 3, ], 
		rmses[2, 1, ], rmses[2, 2, ], rmses[2, 3, ], 
		rmses[3, 1, ], rmses[3, 2, ], rmses[3, 3, ], 
		rmses[4, 1, ], rmses[4, 2, ], rmses[4, 3, ], 
		rmses[5, 1, ], rmses[5, 2, ], rmses[5, 3, ], 
		rmses[6, 1, ], rmses[6, 2, ], rmses[6, 3, ], 
		names = names, las = 2)

names = c()

for (method in METHODS){
	for (c_param in cs){
		names = c(names, paste(method, "c =", c_param))
	}
}

par(mar = c(20,3,3,3))
boxplot(rmses[1, 1, ], rmses[2, 1, ], rmses[3, 1, ], rmses[4, 1, ], rmses[5, 1, ], rmses[6, 1, ], 
		rmses[1, 2, ], rmses[2, 2, ], rmses[3, 2, ], rmses[4, 2, ], rmses[5, 2, ], rmses[6, 2, ],
		rmses[1, 3, ], rmses[2, 3, ], rmses[3, 3, ], rmses[4, 3, ], rmses[5, 3, ], rmses[6, 3, ],
		names = names, las = 2, ylim = c(0, 1.5))

#for each tf, how many genes did it appear in for each method?

gene_by_tf = matrix(0, nrow = MAX_GENE_NUM, ncol = ncol(tf_train))
rownames(gene_by_tf) = 
colnames(gene_by_tf) = colnames(tf_train)
for (g in 1 : MAX_GENE_NUM){
	for (t in 1 : ncol(tf_train)){
		names_of_all_tfs = colnames(tf_train)
		tfs = names(all_results[[g]][['0']][['20']][["important_tfs_at_alpha_simul_max"]])
		cols = which(names_of_all_tfs %in% tfs)
		gene_by_tf[g, cols] = 1
	}
}

colSums(gene_by_tf)
rowSums(gene_by_tf)
sum(colSums(gene_by_tf))