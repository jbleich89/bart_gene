
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
for (j in 1 : length(gene_num_list)){
	gene = gene_num_list[j]
	ybar = mean(gene_train[ , gene ])
	yvec = gene_cv[ , gene]
	control_rmses[j] = sqrt(sum((yvec - ybar) ^ 2) / length(yvec)) 
}

avg_ybar_oormse = mean(control_rmses)

par(mar = c(20,4,3,3))
boxplot(control_rmses, rmses[1, 1, ], rmses[1, 2, ], rmses[1, 3, ], 
		rmses[2, 1, ], rmses[2, 2, ], rmses[2, 3, ], 
		rmses[3, 1, ], rmses[3, 2, ], rmses[3, 3, ], 
		rmses[4, 1, ], rmses[4, 2, ], rmses[4, 3, ], 
		rmses[5, 1, ], rmses[5, 2, ], rmses[5, 3, ], 
		rmses[6, 1, ], rmses[6, 2, ], rmses[6, 3, ], 
		names = names, las = 2, ylim = c(0, 1.8))
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
points(apply(cbind(control_rmses, rmses[1, 1, ], rmses[2, 1, ], rmses[3, 1, ], rmses[4, 1, ], rmses[5, 1, ], rmses[6, 1, ], 
						rmses[1, 2, ], rmses[2, 2, ], rmses[3, 2, ], rmses[4, 2, ], rmses[5, 2, ], rmses[6, 2, ],
						rmses[1, 3, ], rmses[2, 3, ], rmses[3, 3, ], rmses[4, 3, ], rmses[5, 3, ], rmses[6, 3, ]), 2, mean), pch = "-", col = "blue", cex = 5)
abline(h = avg_ybar_oormse, col = "red")