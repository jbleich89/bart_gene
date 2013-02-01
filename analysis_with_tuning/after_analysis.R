tryCatch(library(energy), error = function(e){install.packages("energy")}, finally = library(energy))


#######boxplots
if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/workspace/bart_gene/")
}

source("analysis_with_tuning/simulation_params.R")

setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene")
load("all_results.RData")
load("all_validations.RData")

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






MAX_GENE_NUM = 80

#now we want a boxplot

rmses = array(NA, c(length(cs), length(METHODS), MAX_GENE_NUM))
		
gene_num_list = c(1 : MAX_GENE_NUM)
for (i in 1 : length(cs)){	
	c_param = as.character(cs[i])
	for (j in 1 : length(METHODS)){
		rmses[i, j, ] = sapply(gene_num_list, function(s){all_validations[[s]][[c_param]][["20"]][["0.05"]][[METHODS[j]]]})
	}
}

names = c("Control")
for (c_param in cs){
	for (method in METHODS){
		names = c(names, paste(method, "c =", c_param))
	}
}



##create ybar control
control_rmses = numeric(length(gene_num_list))
for(j in 1 : length(gene_num_list)){
  gene = gene_num_list[j]
  ybar = mean(gene_train[ , gene ])
  yvec = gene_cv[ , gene]
  control_rmses[j] = sqrt(sum((yvec - ybar) ^ 2) / length(yvec)) 
}

avg_ybar_oormse = mean(control_rmses)

par(mar = c(25,8,8,8))
boxplot(control_rmses, rmses[1, 1, ], rmses[1, 2, ], rmses[1, 3, ], 
	rmses[2, 1, ], rmses[2, 2, ], rmses[2, 3, ], 
	rmses[3, 1, ], rmses[3, 2, ], rmses[3, 3, ], 
	rmses[4, 1, ], rmses[4, 2, ], rmses[4, 3, ], 
	rmses[5, 1, ], rmses[5, 2, ], rmses[5, 3, ], 
	rmses[6, 1, ], rmses[6, 2, ], rmses[6, 3, ], 
	names = names, las = 2, ylim = c(0.5, 1.0))
abline(h = avg_ybar_oormse, col = "red")
#add means to the boxlplot
points(apply(cbind(control_rmses, rmses[1, 1, ], rmses[1, 2, ], rmses[1, 3, ], 
	rmses[2, 1, ], rmses[2, 2, ], rmses[2, 3, ], 
	rmses[3, 1, ], rmses[3, 2, ], rmses[3, 3, ], 
	rmses[4, 1, ], rmses[4, 2, ], rmses[4, 3, ], 
	rmses[5, 1, ], rmses[5, 2, ], rmses[5, 3, ], 
	rmses[6, 1, ], rmses[6, 2, ], rmses[6, 3, ]), 2, mean), pch = "-", col = "blue", cex = 5)

names = c("Control")

for (method in METHODS){
	for (c_param in cs){
		names = c(names, paste(method, "c =", c_param))
	}
}

par(mar = c(20,4,3,3))
boxplot(control_rmses, rmses[1, 1, ], rmses[2, 1, ], rmses[3, 1, ], rmses[4, 1, ], rmses[5, 1, ], rmses[6, 1, ], 
	rmses[1, 2, ], rmses[2, 2, ], rmses[3, 2, ], rmses[4, 2, ], rmses[5, 2, ], rmses[6, 2, ],
	rmses[1, 3, ], rmses[2, 3, ], rmses[3, 3, ], rmses[4, 3, ], rmses[5, 3, ], rmses[6, 3, ],
	names = names, las = 2, ylim = c(0, 1.8), 
	main = "Out-of-Sample-RMSEs for all covariate weights\nand all methods to detect importance", 
	ylab = "Out-of-Sample-RMSE")



#for each tf, how many genes did it appear in for each method?

gene_by_tf = matrix(0, nrow = MAX_GENE_NUM, ncol = ncol(tf_train))
rownames(gene_by_tf) = colnames(gene_by_tf) = colnames(tf_train)
for (g in 1 : MAX_GENE_NUM){
	for (t in 1 : ncol(tf_train)){
		names_of_all_tfs = colnames(tf_train)
		tfs = names(all_results[[g]][['0']][['20']][["0.05"]][["important_tfs_at_alpha_pointwise"]])
		cols = which(names_of_all_tfs %in% tfs)
		gene_by_tf[g, cols] = 1
	}
}

colSums(gene_by_tf)
sum(colSums(gene_by_tf))
#% sparsity
sum(rowSums(gene_by_tf) == 0) / nrow(gene_by_tf)

par(mar = c(5,5,5,5))
poissonity_pval = poisson.mtest(rowSums(gene_by_tf), R = 999)$p.value
barplot(table(rowSums(gene_by_tf)), main = paste("# of TF's Affecting all 6,000 Genes\nPoissonity pval =", round(poissonity_pval, 3)), xlab = "# TF's", ylab = "# Genes")

