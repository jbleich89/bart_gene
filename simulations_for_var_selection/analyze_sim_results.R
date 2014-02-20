






setwd("C:/Users/Kapelner/workspace/bart_gene/simulations_for_var_selection/sim_results")

#bar colors for method
model_colors = c(
	rgb(0.1,0.1,0.4),	#BART CV
	rgb(0.3,0.3,1),		#BART ptwise
	rgb(0.5,0.5,1),		#BART simul max
	rgb(0.7,0.7,1),		#BART simul se
	rgb(0.1,0.1,0.2),	#BART CV with GOOD prior
	rgb(0.3,0.3,0.4),	#BART ptwise with GOOD prior
	rgb(0.5,0.5,0.6),	#BART simul max with GOOD prior
	rgb(0.7,0.7,0.8),	#BART simul se with GOOD prior
	rgb(0.3,0.15,0),	#BART CV with BAD prior
	rgb(0.5,0.25,0.0),	#BART ptwise with BAD prior
	rgb(0.7,0.35,0.0),	#BART simul max with BAD prior
	rgb(0.9,0.45,0),		#BART simul se with BAD prior
	rgb(0.9,0.9,0),		#Stepwise B
	rgb(0.9,0.8,0.4),	#Stepwise F
	rgb(1,0.5,0.5),		#Lasso
	rgb(0.1333333, 0.5450980, 0.1333333),		#RF CV
	rgb(0.4,1,0.4),		#RF ptwise
	rgb(0.8,1,0.8)		#RF bonf
)

num_models = length(model_colors)


####FINAL PUB PLOTS
setwd("C:/Users/jbleich/Dropbox/BART_gene/sim_results_new/") ##for justin

CEX_SIZE=2.75
CEX_SIZE_2 = 2
METHOD_NAMES = c("BART-\nBest", "BART-\nLocal", "BART-\nG.Max","BART-\nG.SE", "Step-\nwise", "Lasso\n", "RF-\nCV", "DT\n", "spike-\nslab")

make_lin_plot = function(p, p0, sigsq){
	par(mar = c(6.5,5.5,2,0))
	par(mgp=c(4,5,0))		

	load(file = paste("results_p", p, "_p", p0, "_sigsq", sigsq, ".Rdata", sep = ""))
	obj = eval(as.name(paste("results_p", p, "_p", p0, "_sigsq", sigsq, sep = "")))
	F1s = apply(obj[,,3], 1, mean)[c(1,2,3,4, 13,15,16, 19, 20)]
  
  windows()
	barplot(F1s, 
#		xlab = "Method", 
#		ylab = "F-Score", 
		ylim = c(0, 1),
     		cex.names = CEX_SIZE,
		names = METHOD_NAMES,
		yaxt="n")
#		main = paste("Linear F scores by sigsq and method, p =", p, "and p0 =", p0),
		#col = model_colors[c(1,2,3,4, 13,15,16)])
	
#	axis(1, at = 1:6, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
	title(ylab="F-Score",mgp=c(3.5,1,0), cex.lab = CEX_SIZE)
     axis(2, at= c(0,.2,.4,.6,.8,1), labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),mgp=c(3,1,0),cex.axis=CEX_SIZE_2)
	
}

make_lin_plot(200, 2, 5)
make_lin_plot(200, 2, 20)
make_lin_plot(200, 20, 5)
make_lin_plot(200, 20, 20)


make_nonlin_plot = function(p, sigsq){
	par(mar = c(6.5,5.5,2,0))
	par(mgp=c(4,5,0))	
	
	X = read.csv(paste("complete_var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))
	F1s = X[c(1,2,3,4, 13,15,16), 4]
	barplot(F1s, 
#		xlab = "Method", 
#		ylab = "F-Score",
		yaxt="n", 
		ylim = c(0, 1),
     		cex.names = CEX_SIZE,
		names = METHOD_NAMES,
		yaxt="n",
#		main = paste("Linear F scores by sigsq and method, p =", p, "and p0 =", p0),
			col = model_colors[c(1,2,3,4, 13,15,16)])
	
#	axis(1, at = 1:6, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
	title(ylab="F-Score",mgp=c(3.5,1,0), cex.lab = CEX_SIZE)
     axis(2, at= c(0,.2,.4,.6,.8,1), labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),mgp=c(3,1,0),cex.axis=CEX_SIZE_2)
	
}


#make_nonlin_plot(200, 100)
#make_nonlin_plot(200, 625)
make_nonlin_plot(500, 100)
make_nonlin_plot(500, 625)
make_nonlin_plot(1000, 100)
make_nonlin_plot(1000, 625)

make_lin_prior_plot = function(p, p0, sigsq){
	par(mar = c(6.5,5.5,2,0))
	par(mgp=c(4,5,0))	
	
	X = read.csv(paste("complete_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
	F1s = X[c(1, 5, 9), 4]
	barplot(F1s, 
#		xlab = "method", 
#			ylab = "F-Score", 
			yaxt="n",
			ylim = c(0, 1),
			names = c("Uniform\nPrior", "Correct\nPrior", "Incorrect\nPrior"),
			cex.names = CEX_SIZE,
#		main = paste("Linear F scores by sigsq and method, p =", p, "and p0 =", p0),
			col = model_colors[c(1, 5, 9)])

	title(ylab="F-Score",mgp=c(3.5,1,0), cex.lab = CEX_SIZE)
     	axis(2, at= c(0,.2,.4,.6,.8,1), labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),mgp=c(3,1,0),cex.axis=CEX_SIZE_2)
}

make_lin_prior_plot(200, 2, 5)
make_lin_prior_plot(200, 2, 20)
make_lin_prior_plot(200, 20, 5)
make_lin_prior_plot(200, 20, 20)







setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene/")