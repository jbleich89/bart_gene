
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

##friedman analysis
##friedman analysis
##friedman analysis
ps = c(25, 100, 200, 500, 1000)
sigsqs = c(1, 100, 625)


spacing = 2

graphics.off()
par(mfrow = c(2, 1))
graph_counter = 0
for (p in ps){
	results = array(NA, length(sigsqs) * num_models + spacing * (length(sigsqs) - 1))
	names(results)[seq(from = 7, to = length(sigsqs) * num_models + spacing * (length(sigsqs) - 1), by = length(model_colors))] = sigsqs
	iter = 1
	for (sigsq in sigsqs){
		X = read.csv(paste("complete_var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))
		F1_results = X[, 4]
		F1_results[c(2,3,4,6,7,8,10,11,12,17,18)] = 0
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
sigsqs = c(1, 5, 20)

graphics.off()
par(mfrow = c(2, 2))
graph_counter = 0
master_iter = 0
for (p in ps){
	for (po_prop in po_props){
		p0 = ceiling(p * po_prop)
		
		results = array(NA, length(sigsqs) * num_models + spacing * (length(sigsqs) - 1))
		names(results)[seq(from = 7, to = length(sigsqs) * num_models + spacing * (length(sigsqs) - 1), by = length(model_colors) + spacing)] = sigsqs
		iter = 1
		for (sigsq in sigsqs){
			master_iter = master_iter + 1
			print(master_iter)
			X = read.csv(paste("complete_var_sel_sim_linear_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
			F1_results = X[, 4]
			F1_results[c(2,3,4,6,7,8,10,11,12,17,18)] = 0
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


setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene/writeup_4_23_13")