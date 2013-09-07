##SVD Script 
if (.Platform$OS.type == "windows"){
  directory_where_code_is = "C:\\Users\\jbleich\\workspace\\CGMBART_GPL/bartMachine/R/"
}
setwd(directory_where_code_is)
source("bart_package_builders.R")
source("bart_package_variable_selection.R")
source("bart_package_data_preprocessing.R")
source("bart_package_inits.R")
source("bart_package_summaries.R")


setwd("C:/Users/jbleich/workspace/bart_gene/analysis_with_tuning/")


source("simulation_params_JB.R")

setwd("C:\\Users\\jbleich\\workspace\\CGMBART_GPL/")



N = 250
p = 40
NUM_DATA = 100
NUM_BART = 50
SIGMA = 1
BURN_SIZE = 2000
POST_SIZE = 5000
NTREE = 10

probs_mat_svd = matrix(NA, nrow = NUM_DATA*NUM_BART, ncol = p + 2)
corrs_mat_svd = matrix(NA, nrow = NUM_DATA, ncol = p)
counter = 1

for(i in 1 : NUM_DATA){
  x.train = data.frame(matrix(rnorm(N * p,0,1),nrow = N, ncol = p))
  y = rnorm(N, 0, SIGMA)
  decomp = svd(cbind(y,x.train))
  u = decomp$u
  dim(cor(u))
  uy = u[, 1]
  ux = data.frame(u[, 2 : (p + 1)])
  ##get correlations
  corrs_mat_svd[i, ] = sapply(1 : p, function(s) abs(cor(uy, ux[,s])))
  for(j in 1 : NUM_BART){
    ## do bart stuff
    bart_noise = build_bart_machine(X = ux, y = uy, num_trees = NTREE, num_burn_in = BURN_SIZE, num_iterations_after_burn_in = POST_SIZE, verbose = F, run_in_sample = F)
    avg_props = get_var_props_over_chain(bart_noise)
    destroy_bart_machine(bart_noise)
    
    probs_mat_svd[counter, ] = c(i, j , avg_props)
    
    if(counter %% 25 == 0) print(counter)
    counter = counter + 1
  }
}

save(probs_mat_svd, file = "null_probs_matrix_svd.Rdata")
save(corrs_mat_svd, file = "null_corrs_matrix_svd.Rdata")



#load items
load("C:/Users/Justin/Dropbox/BART_gene/null_corrs_matrix_svd.Rdata")
load("C:/Users/Justin/Dropbox/BART_gene/null_probs_matrix_svd.Rdata")
load("C:/Users/jbleich/Desktop/Dropbox/BART_gene/null_corrs_matrix_svd.Rdata")
load("C:/Users/jbleich/Desktop/Dropbox/BART_gene/null_probs_matrix_svd.Rdata")
setwd("C:/Users/jbleich/Desktop/Dropbox/BART_gene/null_distribution/")
par(mgp=c(1.8,.5,0), mar=c(3.5,3.5,2,1)) 
##across data set variation

var_idx = 3 : 42

sd_by_dataset = sapply(var_idx, function(s) tapply(probs_mat_svd[,s], probs_mat_svd[, 1] ,sd))
hist(apply(sd_by_dataset, 2, mean), breaks = 50, col = "grey", xlab =  "Average SD of Variable Inclusion Proportion", main ="" ) ## avg sd across datasets for each of the genes
hist(c(sd_by_dataset), breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "")

mean_by_dataset = sapply(var_idx, function(s) tapply(probs_mat_svd[,s], probs_mat_svd[, 1] ,mean))
#
hist(mean_by_dataset, breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "")


#overall mean
hist(apply(probs_mat_svd[,var_idx], 2, mean), breaks = 50, col = "grey", xlab =  "Average Variable Inclusion Proportion", main ="" ) 


##within data set variation
sd_by_bart = sapply(var_idx, function(s) tapply(probs_mat_svd[,s], probs_mat_svd[, 2] ,sd))
hist(apply(sd_by_bart, 2, mean), breaks = 50, col = "grey", xlab = "Average SD of Variable Inclusion Proportion", main ="" ) ## avg sd across datasets for each of the genes
hist(c(sd_by_bart), breaks = 50, col = "grey",xlab = "SD of Variable Inclusion Proportion", main = "" )

#Shane comment 3 - dual plot above
##first plot takes sd of numbers across data sets, second across bart runs
##need to stack. 
par(mfrow = c(2,1))
hist(c(sd_by_dataset), breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "", xlim= c(0,.0175)) ##across data set var
hist(c(sd_by_bart), breaks = 50, col = "grey",xlab = "SD of Variable Inclusion Proportion", main = "",xlim= c(0,.0175 )) ##across bart var
#dev.copy2pdf(file = "dual_dataset_bart_svd.pdf", out.type="pdf")
#dev.off()


##shane comment 5:
par(mfrow = c(2,1))
sd_k = apply(probs_mat_svd[,3:42], 1, sd)
hist(c(sd_by_dataset), breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "", xlim= c(0,.025)) ##across data set var     
hist(sd_k, breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "", xlim = c(0,.025))  
#dev.copy2pdf(file = "dual_dataset_var_svd.pdf", out.type="pdf")
#dev.off()

##histogram across data sets for 1 var
#across data sets
sample_var_over_data = tapply(probs_mat_svd[,3], probs_mat_svd[,1], mean)
mean(sample_var_over_data)
hist(sample_var_over_data, breaks = 25, col = "grey", xlab = "Average Variable Inclusion Proportion", main = "")

all_var_over_data = sapply(var_idx, function(s) tapply(probs_mat_svd[,s], probs_mat_svd[, 1] ,mean))  
boxplot(all_var_over_data, ylab = "Variable Inclusion Proprtion", xlab ="Variable") ##Shane comment 1

#dev.copy2pdf(file = "all_vars_across_data_svd.pdf", out.type="pdf")
#dev.off()
#within data sets
sd(sample_var_over_data)
avg_in_bart_run = probs_mat_svd[1 : 50, 3]
hist(avg_in_bart_run)
hist(avg_in_bart_run, breaks = 25, col = "grey", xlab = "Average Variable Inclusion Proportion", main = "")


boxplot(probs_mat_svd[1:50,3:42], ylab = "Variable Inclusion Proprtion", xlab ="Variable") ##shane comment 1 - within single bart over 50 runs     
#dev.copy2pdf(file = "all_vars_across_bart_run_svd.pdf", out.type="pdf")
#dev.off()

#shane comment 2:
s1_data = tapply(probs_mat_svd[,3], probs_mat_svd[,1], mean)
s2_data = tapply(probs_mat_svd[,4], probs_mat_svd[,1], mean)
s1_bart = probs_mat_svd[1:50,3]
s2_bart = probs_mat_svd[1:50,4]
cols = c("red","blue","red","blue")
boxplot(s1_data, s1_bart, s2_data, s2_bart, col = cols,names= c("1","1","2","2"),
        xlab = "Variable", ylab = "Variable Inclusion Proportion")      
#dev.copy2pdf(file = "red_blue_svd.pdf", out.type="pdf")
#dev.off()

##correlation stuff 
avg_by_dataset = sapply(var_idx, function(s) tapply(probs_mat_svd[,s], probs_mat_svd[, 1] ,mean))

dim(avg_by_dataset)

dim(corrs_mat_svd)
hist(corrs_mat_svd, breaks = 40, col = "grey", xlab = "Data Correlations", main = "") ##Shane comment 7.
#dev.copy2pdf(file = "data_corrs_svd.pdf", out.type="pdf")
#dev.off()
var_corrs = sapply(1 : 40, function(s) cor(avg_by_dataset[,s], corrs_mat_svd[,s]))
var_corrs
hist(var_corrs, breaks = 40, col = "grey", xlab = "Correlations", main = "")
#dev.copy2pdf(file = "bart_corrs_svd.pdf", out.type="pdf")
#dev.off()