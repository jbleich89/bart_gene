
setwd("C:/Users/Kapelner/workspace/bart_gene/simulations_for_var_selection/sim_results")

#bar colors for method
model_colors = c(
	rgb(0.1,0.1,0.4),	#BART CV
	rgb(0.3,0.3,1),		#BART ptwise
	rgb(0.5,0.5,1),		#BART simul max
	rgb(0.7,0.7,1),		#BART simul se
	rgb(0.9,0.9,0),		#Stepwise B
	rgb(0.9,0.8,0.4),	#Stepwise F
	rgb(1,0.5,0.5),		#Lasso
	rgb(0.1333333, 0.5450980, 0.1333333),		#RF CV
	rgb(0.4,1,0.4),		#RF ptwise
	rgb(0.8,1,0.8)		#RF bonf
)

num_models = length(model_colors)

##friedman analysis
##friedman analysis
##friedman analysis
ps = c(25, 100, 200, 500, 1000)
sigsqs = c(1, 10, 50, 100)


spacing = 2

graphics.off()
par(mfrow = c(2, 1))
graph_counter = 0
for (p in ps){
	results = array(NA, length(sigsqs) * num_models + spacing * (length(sigsqs) - 1))
	names(results)[seq(from = 6, to = length(sigsqs) * num_models + spacing * (length(sigsqs) - 1), by = 12)] = sigsqs
	iter = 1
	for (sigsq in sigsqs){
		X = read.csv(paste("var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))
		X2 = read.csv(paste("rf_var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))
		F1_results = X[, 4]
		rf_ptwise = X2[8, 4]
		rf_simul = X2[9, 4]		
		rf_cv = X2[10, 4]
		F1_results = c(F1_results[1 : 7], rf_cv, rf_ptwise, rf_simul)
		results[iter : (iter + num_models - 1)] = F1_results
		iter = iter + num_models + spacing
	}
	barplot(results, 
			xlab = "sigsq + method", 
			ylab = "F-score of avg", 
			ylim = c(0, 1),
			main = paste("Friedman F scores by sigsq and method, p =", p),
			col = rep(c(model_colors, rep(NA, spacing)), length(sigsqs)))
	
	graph_counter = graph_counter + 1
	if (graph_counter %% 2 == 0){
		windows()
		par(mfrow = c(2, 1))
	}
}


##linear model analysis
##linear model analysis
##linear model analysis

ps = c(20, 100, 200, 500, 1000)
po_props = c(0.01, 0.05, 0.1, 0.2)
sigsqs = c(0.1, 0.5, 1, 5)

graphics.off()
par(mfrow = c(2, 2))
graph_counter = 0
for (p in ps){
	for (po_prop in po_props){
		p0 = ceiling(p * po_prop)
		
		results = array(NA, length(sigsqs) * num_models + spacing * (length(sigsqs) - 1))
		names(results)[seq(from = 6, to = length(sigsqs) * num_models + spacing * (length(sigsqs) - 1), by = 12)] = sigsqs
		iter = 1
		for (sigsq in sigsqs){
			X = read.csv(paste("var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
			X2 = read.csv(paste("rf_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
#			print(X[,4])
			F1_results = X[, 4]
			rf_ptwise = X2[8, 4]
			rf_simul = X2[9, 4]
			rf_cv = X2[10, 4]
			F1_results = c(F1_results[1:7], rf_cv, rf_ptwise, rf_simul)
			results[iter : (iter + num_models - 1)] = F1_results
			iter = iter + num_models + spacing
		}
		barplot(results, 
				xlab = "sigsq + method", 
				ylab = "F-score of avg", 
				ylim = c(0, 1),
				main = paste("Linear F scores by sigsq and method, p =", p, "and p0 =", p0),
				col = rep(c(model_colors, rep(NA, spacing)), length(sigsqs)))
		
		graph_counter = graph_counter + 1
		if (graph_counter %% 4 == 0){
			windows()
			par(mfrow = c(2, 2))
		}		
	}

}


setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene/writeup_4_3_13")