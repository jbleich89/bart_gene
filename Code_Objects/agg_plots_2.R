require(BayesTree)
setwd("~/Documents/Research/Genomics")
source("bart_fns.R")
source("bart_fns2.R")
source("agg_plots_fns.R")

#load("bart301_v3.R")
#load("bart501_v3.R")
#load("bart1001_v3.R")
#load("bart2001_v3.R")
#load("bart3001_v3.R")

set.seed(20)

##TF settings
tf.size1=c(10,15,20,25)
tf.size2=c(10,20,30,40) ##script is set to work with this right now
tfsize.vec=tf.size2
n.tree.vec=c(5,10,15,20)
ntree.vec=n.tree.vec
percent=0
fn_name="median"
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12)

##Code for having TFs on X-axis (hold tf const)
##Matrix and plot

par(mfrow=c(1,1))

tf.const.mat=matrix(nrow=5,ncol=(length(tf.size2)+1))
dim(tf.const.mat)

tf.const.mat[1,]=c(30,tf_axis_aggregate(output.30.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tf.const.mat[2,]=c(50,tf_axis_aggregate(output.50.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tf.const.mat[3,]=c(100,tf_axis_aggregate(output.100.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tf.const.mat[4,]=c(200,tf_axis_aggregate(output.200.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tf.const.mat[5,]=c(300,tf_axis_aggregate(output.300.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tf.const.mat
bart_agg_plot_tf_axis(tf.const.mat,tf.size2,factor.vec2)

##Code for Having Trees on X-axis (hold tree constant)
##Create matrix
tree.const.mat=matrix(nrow=5,ncol=(length(n.tree.vec)+1))
dim(tree.const.mat)
tree.const.mat
tree.const.mat[1,]=c(30,tree_axis_aggregate(output.30.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tree.const.mat[2,]=c(50,tree_axis_aggregate(output.50.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tree.const.mat[3,]=c(100,tree_axis_aggregate(output.100.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tree.const.mat[4,]=c(200,tree_axis_aggregate(output.200.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tree.const.mat[5,]=c(300,tree_axis_aggregate(output.300.1,tf.size2,ntree.vec=n.tree.vec,0,fn_name))
tree.const.mat
bart_agg_plot_tree_axis(tree.const.mat,n.tree.vec,factor.vec2)




##Junk
calc_across_tfs_for_plot(output.30.1,tf.size2,ntree.vec=n.tree.vec,0,"median")
bart_plot_both(output=output.200.1,ntree.vec=n.tree.vec,tfsize.vec=tf.size2,factor.vec=factor.vec2,num.obs=30,percent=0)