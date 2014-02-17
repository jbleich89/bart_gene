

#run on the GRID ONLY!!
setwd("../CGMBART_GPL")

MAX_GENE_NUM = 6026

simulation_names = c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso-CV", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full", "Rob-Best", "dynaTree", "spikeslab")
all_rmse_results = matrix(NA, nrow = MAX_GENE_NUM, ncol = length(simulation_names))
colnames(all_rmse_results) = simulation_names
rownames(all_rmse_results) = 1 : MAX_GENE_NUM
all_num_var_results = matrix(NA, nrow = MAX_GENE_NUM, ncol = length(simulation_names))
colnames(all_num_var_results) = simulation_names
rownames(all_num_var_results) = 1 : MAX_GENE_NUM


for (g in 1 : MAX_GENE_NUM){
	if (g %in% c(4124, 4127, 4130, 4133, 4135, 4136, 4137, 4138, 4139, 6026)){
		next
	}
	load(file = paste("rmse_results_", g, ".RData", sep = ""))
	load(file = paste("num_var_results_", g, ".RData", sep = ""))
	all_rmse_results[g, ] = rmse_results
	rownames(all_rmse_results)[g] = rownames(rmse_results)
	all_num_var_results[g, ] = num_var_results
	rownames(all_num_var_results)[g] = rownames(num_var_results)[1]
}

head(all_rmse_results)
head(all_num_var_results)
tail(all_rmse_results)
tail(all_num_var_results)

save(all_rmse_results, file = paste("all_rmse_results.RData", sep = ""))
save(all_num_var_results, file = paste("all_num_var_results.RData", sep = ""))

###Get DT and Spike

setwd("C:/Users/jbleich/Dropbox/BART_gene/dt_spike_results_raw")

MAX_GENE_NUM = 6026

simulation_names = c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso-CV", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full", "Rob-Best", "dynaTree", "spikeslab")
all_rmse_results_dt_spike = matrix(NA, nrow = MAX_GENE_NUM, ncol = length(simulation_names))
colnames(all_rmse_results_dt_spike) = simulation_names
rownames(all_rmse_results_dt_spike) = 1 : MAX_GENE_NUM
all_num_var_results_dt_spike = matrix(NA, nrow = MAX_GENE_NUM, ncol = length(simulation_names))
colnames(all_num_var_results_dt_spike) = simulation_names
rownames(all_num_var_results_dt_spike) = 1 : MAX_GENE_NUM


for (g in 1 : MAX_GENE_NUM){
  if (g %in% c(4124, 4127, 4130, 4133, 4135, 4136, 4137, 4138, 4139, 6026)){
    next
  }
  rmse_results = matrix(rep(NA,12), ncol = 13)
  num_var_results = NA
  tryCatch(load(file = paste("rmse_results_dt_spike_", g, ".RData", sep = "")), error = function(e) {return(1)})
  tryCatch(load(file = paste("num_var_results_dt_spike_", g, ".RData", sep = "")), error = function(e) {print(g)})
  if(is.na(rmse_results[,12]) ) next
  all_rmse_results_dt_spike[g, ] = rmse_results
  rownames(all_rmse_results_dt_spike)[g] = rownames(rmse_results)
  all_num_var_results_dt_spike[g, ] = num_var_results
  rownames(all_num_var_results_dt_spike)[g] = rownames(num_var_results)[1]
}


median(all_rmse_results_dt_spike[,12], na.rm = T)


head(all_rmse_results_dt_spike)
head(all_num_var_results_dt_spike)
boxplot(all_num_var_results_dt_spike[,12:13])
tail(all_rmse_results)
tail(all_num_var_results)


redos = unique(c(which(is.na(all_rmse_results_dt_spike[,12])),
which(all_rmse_results_dt_spike[,12] > 10)))

max(all_rmse_results_dt_spike[,12], na.rm = T)

paste(redos, collapse = ",")


new_all_rmse_results_wanted = cbind(all_rmse_results_wanted, all_rmse_results_dt_spike[,12:13])
head(new_all_rmse_wanted)
boxplot(new_all_rmse_results_wanted, na.rm = T)
save(new_all_rmse_results_wanted, file = "new_all_rmse_results_wanted.Rdata")

new_all_num_var_results_no_na_wanted = cbind(all_num_var_results_no_na_wanted, all_num_var_results_dt_spike[,12 :13])
boxplot(new_all_num_var_results_no_na_wanted, na.rm = T)
save(new_all_num_var_results_no_na_wanted, file = "new_all_num_var_results_no_na_wanted.Rdata")

#########load results into R at this point ON THE LOCAL MACHINE
##work
#load("C:/Users/jbleich/Dropbox/BART_gene/all_rmse_results_wanted.RData")
 #load("C:/Users/jbleich/Dropbox/BART_gene/all_num_var_results.RData")
 #load("C:/Users/jbleich//Dropbox/BART_gene/all_num_var_results_no_na_wanted.RData")
# load("C:/Users/jbleich/Desktop/Dropbox/BART_gene/all_rmse_per_num_var_wanted.RData")
# #home
# load("C:/Users/Justin/Dropbox/BART_gene/all_rmse_per_num_var_wanted.RData")
# load("C:/Users/Justin/Dropbox/BART_gene/all_num_var_results_no_na_wanted.RData")
# load("C:/Users/Justin/Dropbox/BART_gene/all_rmse_results_wanted.RData")
# load("C:/Users/Justin/Dropbox/BART_gene/all_rmse_results.RData")
# load("C:/Users/Justin/Dropbox/BART_gene/all_num_var_results.RData")

load("C:/Users/jbleich/Dropbox/BART_gene/new_all_rmse_results_wanted.Rdata")
load("C:/Users/jbleich/Dropbox/BART_gene/new_all_num_var_results_no_na_wanted.Rdata")


##Change some colnames
colnames(new_all_num_var_results_no_na_wanted)[8]="BART-G.Max"
colnames(new_all_rmse_per_num_var_wanted)[8]="BART-G.Max"
colnames(new_all_rmse_results_wanted)[8]="BART-G.Max"

colnames(new_all_num_var_results_no_na_wanted)[10]="DT"
colnames(new_all_rmse_results_wanted)[10]="DT"

colnames(new_all_num_var_results_no_na_wanted)[11]="spike-slab"
colnames(new_all_rmse_results_wanted)[11]="spike-slab"

rm = which(new_all_num_var_results_no_na_wanted[,1]==39) ##why? these are bad genes we didn't get results for
all_num_var_results_no_na_wanted = new_all_num_var_results_no_na_wanted[-rm,]
#all_rmse_per_num_var_wanted = all_rmse_per_num_var_wanted[-rm,]
all_rmse_results_wanted = new_all_rmse_results_wanted[-rm,]

#oos rmses for all methods


#all_rmse_results_wanted = all_rmse_results[, c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full")]
#save(all_rmse_results_wanted, file = "all_rmse_results_wanted.RData")

##rm OLS BART BEST
all_rmse_results_wanted = all_rmse_results_wanted[, c("Null", "OLS", "Stepwise", "Lasso", "RF", "BART-Best", "BART-G.Max", "BART-Full", "DT", "spike-slab")]

par(mgp=c(1.8,.5,0),mar=c(3,3,2,1))
boxplot(all_rmse_results_wanted, ylim = c(0, 1.1), ylab = "out-of-sample RMSE", outline = F)
points(apply(all_rmse_results_wanted, 2, mean, na.rm = T), pch = "-", col = "blue", cex = 5)
#abline(a = mean(all_rmse_results_wanted[, 1], na.rm = TRUE), b = 0, col = "red")


##NUM TFs selected###################

#all_num_var_results_no_na = t(sapply(1 : nrow(all_num_var_results), function(i){ifelse(is.na(all_num_var_results[i, ]), 39, all_num_var_results[i, ])}))
#all_num_var_results_no_na_wanted = all_num_var_results_no_na[, c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso", "RF", "BART-Best", "BART-G.Max", "BART-Full")]
#save(all_num_var_results_no_na_wanted, file = "all_num_var_results_no_na_wanted.RData")
#save(all_num_var_results_no_na_wanted, file = "all_num_var_results_no_na_wanted.RData")

##rm ols bart best
all_num_var_results_no_na_wanted = all_num_var_results_no_na_wanted[, c("Null", "OLS", "Stepwise", "Lasso", "RF", "BART-Best", "BART-G.Max", "BART-Full", "DT", "spike-slab")]

par(mgp=c(1.8,.5,0),mar=c(3,3,2,1))
boxplot(all_num_var_results_no_na_wanted, ylab = "Number of TF's Selected", outline = F)
points(apply(all_num_var_results_no_na_wanted, 2, mean), pch = "-", col = "blue", cex = 5)


##Reduction per var 
all_rmse_minus_null_per_num_var = (matrix(rep(all_rmse_results_wanted[, 1], ncol(all_rmse_results_wanted)), ncol = ncol(all_rmse_results_wanted)) - all_rmse_results_wanted) / all_num_var_results_no_na_wanted

for (i in 1 : nrow(all_rmse_results_wanted)){
	all_rmse_minus_null_per_num_var[i, ] = ifelse(is.infinite(all_rmse_minus_null_per_num_var[i, ]), NA, all_rmse_minus_null_per_num_var[i, ])
}
#

all_rmse_minus_null_per_num_var_wanted = all_rmse_minus_null_per_num_var[, c("OLS", "Stepwise", "Lasso", "RF", "BART-Best", "BART-G.Max", "BART-Full", "DT", "spike-slab")]

par(mgp=c(1.8,.5,0),mar=c(3,3,2,1))
boxplot(all_rmse_minus_null_per_num_var_wanted, 
		ylab = "RMSE Reduction per TF", 
		ylim = c(-.05, 0.2), outline = F)
points(apply(all_rmse_minus_null_per_num_var_wanted, 2, mean, na.rm = TRUE), pch = "-", col = "blue", cex = 5)
abline(a = 0, b = 0, col = "red")

counts = apply(all_num_var_results_no_na_wanted, 2, function(s) sum(s>0))
counts


##############################################################3
##### rmse per num vars- old

all_rmse_per_num_var = all_rmse_results / all_num_var_results_no_na
for (i in 1 : nrow(all_rmse_results)){
  all_rmse_per_num_var[i, ] = ifelse(is.infinite(all_rmse_per_num_var[i, ]), NA, all_rmse_per_num_var[i, ])
}


all_rmse_per_num_var_wanted = all_rmse_per_num_var[, c("Null", "OLS", "OLS-BART-Best", "Stepwise", "Lasso", "RF", "BART-Best", "BART-S.Max", "BART-Full")]
#save(all_rmse_per_num_var_wanted, file = "all_rmse_per_num_var_wanted.RData")
par(mgp=c(1.8,.5,0),mar=c(3,3,2,1))
boxplot(all_rmse_per_num_var_wanted, ylab = "RMSE / variable", main = "Out-of-Sample RMSE per TF by Method\nConditional on Finding At Least One Variable", ylim = c(0, 0.5))
points(apply(all_rmse_per_num_var_wanted, 2, mean, na.rm = TRUE), pch = "-", col = "blue", cex = 5)
abline(a = mean(all_rmse_per_num_var_wanted[, 1], na.rm = TRUE), b = 0, col = "red")

