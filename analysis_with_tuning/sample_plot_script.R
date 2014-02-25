library(bartMachine)

setwd("C:/Users/jbleich/workspace/bart_gene/analysis_with_tuning/")

source("simulation_params_JB.R")

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



#load items
load("C:/Users/Justin/Dropbox/BART_gene/null_corrs_matrix.Rdata")
load("C:/Users/Justin/Dropbox/BART_gene/null_probs_matrix.Rdata")
load("C:/Users/jbleich/Dropbox/BART_gene/null_corrs_matrix.Rdata")
load("C:/Users/jbleich/Dropbox/BART_gene/null_probs_matrix.Rdata")
setwd("C:/Users/jbleich/Dropbox/BART_gene/null_distribution/")
par(mgp=c(1.8,.5,0), mar=c(3.5,3.5,2,1)) 
##across data set variation

var_idx = 3 : 42
     
sd_by_dataset = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,sd))
hist(apply(sd_by_dataset, 2, mean), breaks = 50, col = "grey", xlab =  "Average SD of Variable Inclusion Proportion", main ="" ) ## avg sd across datasets for each of the genes
hist(c(sd_by_dataset), breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "")

mean_by_dataset = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,mean))
#
hist(mean_by_dataset, breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "")


#overall mean
hist(apply(probs_mat[,var_idx], 2, mean), breaks = 50, col = "grey", xlab =  "Average Variable Inclusion Proportion", main ="" ) 


##within data set variation
sd_by_bart = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 2] ,sd))
hist(apply(sd_by_bart, 2, mean), breaks = 50, col = "grey", xlab = "Average SD of Variable Inclusion Proportion", main ="" ) ## avg sd across datasets for each of the genes
hist(c(sd_by_bart), breaks = 50, col = "grey",xlab = "SD of Variable Inclusion Proportion", main = "" )

#Shane comment 3 - dual plot above
##first plot takes sd of numbers across data sets, second across bart runs
##need to stack. 
par(mfrow = c(2,1))
hist(c(sd_by_dataset), breaks = 200, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "", xlim= c(0,.0175)) ##across data set var
hist(c(sd_by_bart), breaks = 50, col = "grey",xlab = "SD of Variable Inclusion Proportion", main = "",xlim= c(0,.0175 )) ##across bart var
dev.copy2pdf(file = "dual_dataset_bart.pdf", out.type="pdf")
dev.off()


##shane comment 5:
par(mfrow = c(2,1))
sd_k = apply(probs_mat[,3:42], 1, sd)
hist(c(sd_by_dataset), breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "", xlim= c(0,.025)) ##across data set var     
hist(sd_k, breaks = 50, col = "grey", xlab = "SD of Variable Inclusion Proportion", main = "", xlim = c(0,.025))  
#dev.copy2pdf(file = "dual_dataset_var.pdf", out.type="pdf")
#dev.off()

##histogram across data sets for 1 var
#across data sets
sample_var_over_data = tapply(probs_mat[,3], probs_mat[,1], mean)
mean(sample_var_over_data)
hist(sample_var_over_data, breaks = 25, col = "grey", xlab = "Average Variable Inclusion Proportion", main = "")

all_var_over_data = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,mean))  
boxplot(all_var_over_data, ylab = "Variable Inclusion Proprtion", xlab ="Variable") ##Shane comment 1

#dev.copy2pdf(file = "all_vars_across_data.pdf", out.type="pdf")
#dev.off()
     #within data sets
sd(sample_var_over_data)
avg_in_bart_run = probs_mat[1 : 50, 3]
hist(avg_in_bart_run)
hist(avg_in_bart_run, breaks = 25, col = "grey", xlab = "Average Variable Inclusion Proportion", main = "")


boxplot(probs_mat[1:50,3:42], ylab = "Variable Inclusion Proprtion", xlab ="Variable") ##shane comment 1 - within single bart over 50 runs     
#dev.copy2pdf(file = "all_vars_across_bart_run.pdf", out.type="pdf")
#dev.off()

#shane comment 2:
s1_data = tapply(probs_mat[,3], probs_mat[,1], mean)
s2_data = tapply(probs_mat[,4], probs_mat[,1], mean)
s1_bart = probs_mat[1:50,3]
s2_bart = probs_mat[1:50,4]
cols = c("red","blue","red","blue")
boxplot(s1_data, s1_bart, s2_data, s2_bart, col = cols,names= c("1","1","2","2"),
        xlab = "Variable", ylab = "Variable Inclusion Proportion")      
dev.copy2pdf(file = "red_blue.pdf", out.type="pdf")
dev.off()

##correlation stuff 
avg_by_dataset = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,mean))

dim(avg_by_dataset)

dim(corrs_mat)
hist(corrs_mat, breaks = 40, col = "grey", xlab = "Data Correlations", main = "") ##Shane comment 7.
#dev.copy2pdf(file = "data_corrs.pdf", out.type="pdf")
#dev.off()
var_corrs = sapply(1 : 40, function(s) cor(avg_by_dataset[,s], corrs_mat[,s]))
var_corrs
hist(var_corrs, breaks = 40, col = "grey", xlab = "Correlations", main = "")

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


##ANOVA
anova_mat = matrix(nrow = 1, ncol = 3)
for(i in 1:nrow(probs_mat)){
  temp = matrix(rep(probs_mat[i,1:2], each = 40), ncol =2)
  anova_mat = rbind(anova_mat,cbind(temp,probs_mat[i,3:42]))
  if(i %% 50 == 0) print(i)
}

anova_mat = anova_mat[-1,]     
dim(anova_mat)
anova_data = data.frame(anova_mat)
anova_data$X1 = as.factor(anova_data$X1)
anova_data$X2 = as.factor(anova_data$X2)
colnames(anova_data) = c("dataset","bart_run","inclusion_prop")
str(anova_data)
anova_mod = lm(inclusion_prop ~ dataset + bart_run,data = anova_data)
summary(anova_mod)
write.table(anova_data,file="anova_data.txt")



##Nested Anova calcs
#var i_k
p_i_dot_k = sapply(var_idx, function(s) tapply(probs_mat[,s], probs_mat[, 1] ,mean))
dim(p_i_dot_k)
var_i_dot_k = matrix(nrow = NUM_DATA, ncol = p)
for(i in 1 : NUM_DATA){
  
  temp = probs_mat[which(probs_mat[ ,1] == i), var_idx] ##just get what we need
  var_i_dot_k[i,] = sapply(1:p, function(s){ sum((temp[,s] - p_i_dot_k[i ,(s)])^2)/50})
}

##var_k
p_dot_dot_k = apply(probs_mat[ ,var_idx], 2, mean)
var_k = sapply(1:p, function(s) {sum((p_i_dot_k[ ,s] - p_dot_dot_k[s])^2)/100} )
var_k


##var
grand_mean = mean(apply(probs_mat[,var_idx], 2, mean))
grand_mean
col_means = apply(probs_mat[, var_idx],2, mean)
var = sum((col_means - grand_mean)^2)/p
var
     


par(mgp=c(3.5,.5,0), mar=c(5,5,2,1)) 
boxplot(sqrt(var_i_dot_k), ylim = c(0,.017), ylab = "Standard Deviation", xlab = "Predictor", las = 1, outline = F)
points(sqrt(var_k), pch = 16, col = "blue", cex = 2.5)
abline(h = sqrt(var), col = "red", lty = 2, lwd = 3)
dev.copy2pdf(file="anova_plot.pdf", out.type = "pdf")
dev.off()

##scratch
samp = probs_mat[1:50,3:42]
samp_mean = apply(samp,2,mean)
sapply(1:ncol(samp), function(s){ sum((samp[,s]-samp_mean[s])^2)/50 })
