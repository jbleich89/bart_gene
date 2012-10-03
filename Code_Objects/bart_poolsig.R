require(BayesTree)
 setwd("~/Research/Genomics")
source("bart_fns.R")
g=read.table("expression.genes.txt",header=T)
t=read.table("expression.tfs.39.txt",header=T)
dim(g)
dim(t)
head(g)
head(t)
##Set-up
set.seed(20)

tf.size=c(10,15,20,25)
n.gene=75
n.obs1=30
n.obs2=50

mean.tf=0 #1.4862e-05 ##generate the X matrix which will be fixed- sl
sd.tf=1 #.4778817-not using just to keep larger numbers at play

tf.mat.30=sapply(rep(n.obs1*n.gene,max(tf.size)),rnorm,mean=mean.tf,sd=sd.tf) #gives full matrix
tf.mat.50=sapply(rep(n.obs2*n.gene,max(tf.size)),rnorm,mean=mean.tf,sd=sd.tf)
##just reduce for consistency if you only want to experiment with less tfs


tf.beta.1=c(1,rep(0,times=max(tf.size-1)))
tf.beta.2=c(2,rep(0,times=max(tf.size-1)))