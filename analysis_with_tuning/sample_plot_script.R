
if (.Platform$OS.type == "windows"){
  directory_where_code_is = "C:\\Users\\jbleich\\workspace\\CGMBART_GPL"
}
setwd(directory_where_code_is)
source("r_scripts/bart_package_inits.R")
source("r_scripts/bart_package_data_preprocessing.R")
source("r_scripts/bart_package_builders.R")
source("r_scripts/bart_package_plots.R") ##altered to trash label at the top
source("r_scripts/bart_package_variable_selection.R")
source("r_scripts/bart_package_f_tests.R")

setwd("C:/Users/kapelner/workspace/bart_gene/analysis_with_tuning/")
      

source("simulation_params_JB.R")

setwd(directory_where_code_is)

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