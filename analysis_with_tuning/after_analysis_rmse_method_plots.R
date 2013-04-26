


gene_names2 = gene_names[-c(4124, 4127, 4130, 4133, 4135, 4136, 4137, 4138, 4139)]
rmses = array(NA, c(length(cs), length(METHODS), length(gene_names2)))


for (i in 1 : length(cs)){	
	c_param = as.character(cs[i])
	for (j in 1 : length(METHODS)){
		cat("c_param", c_param, "method", METHODS[j], "\n")
		rmse_vec = sapply(gene_names2, function(s){
#					cat("s", s, "c_param", c_param, "method", METHODS[j], "\n")
					all_validations[[s]][[c_param]][["20"]][["0.05"]][[METHODS[j]]]})
		if (class(rmse_vec) == "numeric"){
			rmses[i, j, ] = rmse_vec
		} 
		else {
			rmses[i, j, ] = as.matrix(rmse_vec)
		}
		
	}
}

names = c("Control")
for (c_param in cs){
	for (method in METHODS){
		names = c(names, paste(method, "c =", c_param))
	}
}

##create ybar control
control_rmses = numeric(length(gene_names2))
for (j in 1 : length(gene_names2)){
	gene = gene_names2[j]
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