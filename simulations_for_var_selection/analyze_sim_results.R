
setwd("C:/Users/Kapelner/workspace/bart_gene/simulations_for_var_selection/sim_results")

#bar colors for method
colors = c(
	rgb(0.3,0.3,1),
	rgb(0.5,0.5,1),
	rgb(0.7,0.7,1),
	rgb(0.4,1,0.4),
	rgb(0.7,1,0.7),
	rgb(1,0.5,0.5),
	rgb(0,0,0),
	rgb(0,0,0)	
)

##friedman analysis
##friedman analysis
##friedman analysis
ps = c(25, 100, 200, 500, 1000, 5000)
sigsqs = c(0.1, 0.5, 1, 5)



graphics.off()
par(mfrow = c(2, 1))
graph_counter = 0
for (p in ps){
	results = array(NA, length(sigsqs) * 6 + 7)
	names(results)[seq(3, length(sigsqs) * 6 + 7, 8)] = sigsqs
	iter = 1
	for (sigsq in sigsqs){
		X = read.csv(paste("var_sel_sim_friedman_p", p, "_sigsq_", sigsq, ".csv", sep = ""))
		print(X[,4])
		results[iter : (iter + 5)] = X[, 4]
		iter = iter + 8
	}
	barplot(results, 
			xlab = "sigsq + method", 
			ylab = "avg F-score", 
			ylim = c(0, 1),
			main = paste("F scores by sigsq and method, p =", p),
			col = rep(colors, length(sigsqs)))
	
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
		
		results = array(NA, length(sigsqs) * 6 + 7)
		names(results)[seq(3, length(sigsqs) * 6 + 7, 8)] = sigsqs
		iter = 1
		for (sigsq in sigsqs){
			X = read.csv(paste("var_sel_sim_p", p, "_p0_", p0, "_sigsq_", sigsq, ".csv", sep = ""))
			print(X[,4])
			results[iter : (iter + 5)] = X[, 4]
			iter = iter + 8
		}
		barplot(results, 
				xlab = "sigsq + method", 
				ylab = "avg F-score", 
				ylim = c(0, 1),
				main = paste("F scores by sigsq and method, p =", p, "and p0 =", p0),
				col = rep(colors, length(sigsqs)))
		
		graph_counter = graph_counter + 1
		if (graph_counter %% 4 == 0){
			windows()
			par(mfrow = c(2, 2))
		}		
	}

}