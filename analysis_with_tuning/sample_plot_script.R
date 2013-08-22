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


##For sample figure

y = gene_train[,1]
x = data.frame(tf_train)


bart_machine = build_bart_machine(X=x, y=y)

investigate_var_importance(bart_machine, num_replicates_for_avg = 1, num_var_plot = 15)
counts = get_var_counts_over_chain(bart_machine)

sel = var_selection_by_permute_response_three_methods(bart_machine, num_reps_for_avg = 5, num_var_plot = 20)


##NOISE PLOT
##
out = list()
par(mgp = c(1.8,.5,0),mar = c(3,3,2,1))
par(mfrow = c(1,2))
nsim = 50
ymin = 0
N = 250
p = 40

nrows = N
ncols = p
x.train = data.frame(matrix(rnorm(nrows*ncols,0,1),nrow=nrows,ncol=ncols))
sig_list = c(.1,.5,1,2,4,8)
sigma = sig_list[3]
sigma
n_iter = 12000

y = rnorm(nrows,0,sigma)
#bart.noise=bart(x.train=x.train,y.train=y,ntree=10,nskip = 1,ndpost = n_iter,keepevery=1)
bart.noise = build_bart_machine(X = x.train, y = y, num_trees = 10, num_burn_in = 1, num_iterations_after_burn_in = n_iter)
counts = get_var_counts_over_chain(bart.noise)

cums = apply(counts,2,cumsum)
cum_props = cums/apply(cums,1,sum)
ymax = max(cum_props[2000:n_iter,])
ymax
rain = rainbow(p)

plot(1:n_iter,cum_props[,1],type="l",col=rain[1],ylim=c(ymin,ymax),
     xlab="Gibbs Sample",ylab="Variable Inclusion Proportion")
sapply(2:p, function(s) points(1:n_iter,cum_props[,s],type="l",col=rain[s]))
abline(h = .025, col = "black", lwd = 3, lty = 2)
##follow max
maxs=apply(cum_props,1,max)
mins=apply(cum_props,1,min)
plot(1:n_iter,maxs,ylim=c(ymin,ymax),type="l",col="red",lwd=3,ylab="Variable Inclusion Proportion",xlab="Gibbs Sample")
points(1:n_iter,mins,type="l",col="blue",lwd=3)
abline(h = .025, col = "black", lwd = 3, lty = 2)



####NOISE SIMULATIONS#####################
N = 250
p = 40
NUM_DATA = 100
NUM_BART = 50
SIGMA = 1
BURN_SIZE = 2000
POST_SIZE = 5000
NTREE = 10

probs_mat = matrix(NA, nrow = NUM_DATA*NUM_BART, ncol = p + 2)
corrs_mat = matrix(NA, nrow = NUM_DATA, ncol = p)
counter = 1

for(i in 1 : NUM_DATA){
  x.train = data.frame(matrix(rnorm(N * p,0,1),nrow = N, ncol = p))
  y = rnorm(N, 0, SIGMA)
  ##get correlations
  corrs_mat[i, ] = sapply(1 : p, function(s) abs(cor(y, x.train[,s])))
  for(j in 1 : NUM_BART){
    ## do bart stuff
    bart_noise = build_bart_machine(X = x.train, y = y, num_trees = NTREE, num_burn_in = BURN_SIZE, num_iterations_after_burn_in = POST_SIZE, verbose = F, run_in_sample = F)
    avg_props = get_var_props_over_chain(bart_noise)
    destroy_bart_machine(bart_noise)
    
    probs_mat[counter, ] = c(i, j , avg_props)
    
    if(counter %% 25 == 0) print(counter)
    counter = counter + 1
  }
}

save(probs_mat, file = "null_probs_matrix.Rdata")
save(corrs_mat, file = "null_corrs_matrix.Rdata")

var_idx = 3 : 42

#load items
load("C:/Users/jbleich/Desktop/Dropbox/BART_gene/null_corrs_matrix.Rdata")
load("C:/Users/jbleich/Desktop/Dropbox/BART_gene/null_probs_matrix.Rdata")
par(mgp=c(1.8,.5,0), mar=c(3.5,3.5,2,1)) 
##within data set variation

sd_by_dataset = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,sd))
hist(apply(sd_by_dataset, 2, mean), breaks = 50, col = "grey", xlab =  "Average SD of Variable Inclusion Proportion", main ="" ) ## avg sd across datasets for each of the genes
hist(c(sd_by_dataset), breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "")

##across data set variation
sd_by_bart = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 2] ,sd))
hist(apply(sd_by_bart, 2, mean), breaks = 50, col = "grey", xlab = "Average SD of Variable Inclusion Proportion", main ="" ) ## avg sd across datasets for each of the genes
hist(c(sd_by_bart), breaks = 50, col = "grey",xlab = "SD of Variable Inclusion Proportion", main = "" )


##histogram across data sets for 1 var
#across data sets
sample_var_over_data = tapply(probs_mat[,3], probs_mat[,1], mean)
mean(sample_var_over_data)
hist(sample_var_over_data, breaks = 25, col = "grey", xlab = "Average Variable Inclusion Proportion", main = "")
#within data sets
sd(sample_var_over_data)
avg_in_bart_run = probs_mat[1 : 50, 3]
hist(avg_in_bart_run)
hist(avg_in_bart_run, breaks = 25, col = "grey", xlab = "Average Variable Inclusion Proportion", main = "")

##correlation stuff 
avg_by_dataset = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,mean))

dim(avg_by_dataset)

dim(corrs_mat)
var_corrs = sapply(1 : 40, function(s) cor(avg_by_dataset[,s], corrs_mat[,s]))
var_corrs
hist(var_corrs, breaks = 20, col = "grey")

###SVD for orthogonal stuff
x.train = data.frame(matrix(rnorm(N * p,0,1),nrow = N, ncol = p))
y = rnorm(N, 0, SIGMA)
decomp = svd(cbind(y,x.train))
u = decomp$u
dim(cor(u))
uy = u[, 1]
ux = data.frame(u[, 2 : (p + 1)])

nsim_svd = 1
svd_mat = matrix(nrow = nsim_svd, ncol = p)

for(j in 1 : nsim_svd){
  bart_svd = build_bart_machine(X = ux, y = uy, num_trees = NTREE, num_burn_in = 1, num_iterations_after_burn_in = POST_SIZE, verbose = F, run_in_sample = F)
#  sum_var_counts = get_var_counts_over_chain(bart_svd)
  avg_var_counts_by_col = colSums(var_counts)
  total_count = sum(sum_var_counts_by_col)
  svd_mat[j, ] = get_var_props_over_chain(bart_svd)
  #destroy_bart_machine(bart_svd)
}


counts = get_var_counts_over_chain(bart_svd)

cums = apply(counts,2,cumsum)
cum_props = cums/apply(cums,1,sum)
dim(cum_props)
ymax = max(cum_props[100 : POST_SIZE,])
ymin = min(cum_props[10 : POST_SIZE,])
ymax
rain = rainbow(p)

plot(1 : (POST_SIZE), cum_props[,1], type = "l", col = rain[1], ylim = c(ymin , ymax),
     xlab = "Gibbs Sample", ylab = "Variable Inclusion Proportion")
sapply(2 : p, function(s) points(1 : ( POST_SIZE) ,cum_props[,s], type = "l", col = rain[s]))
abline(h = .025, col = "black", lwd = 3, lty = 2)

#exp_count = rep(total_count/p, times = p)
#chisq.test(sum_var_counts_by_col, p = rep(.025, p))
#total_count

colMeans(svd_mat)
hist(colMeans(svd_mat), col = "grey", main = "", breaks = 25, xlab = "Variable Inclusion Proportion")

###could be non-linear correlation






