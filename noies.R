setwd("~/Research_Genomics/bart_gene") ##work
setwd("~/Documents/Research/Genomics/bart_gene/") ##home
source("real_functions.R")
library(BayesTree)

##work
priors=read.table("~/Research_Genomics/CHIP.priorprobs.39.txt",header=T) ##work
gene.exp=read.table("~/Research_Genomics/expression.genes.txt",header=T) ##work
tf.exp=read.table("~/Research_Genomics/expression.tfs.39.txt",header=T) ##work
##Home
priors=read.table("~/Documents/Research/Genomics/Real_Data/CHIP.priorprobs.39.txt",header=T) ##home
gene.exp=read.table("~/Documents/Research/Genomics/Real_Data/expression.genes.txt",header=T) ##home
tf.exp=read.table("~/Documents/Research/Genomics/Real_Data/expression.tfs.39.txt",header=T) ##home

#setup
priorMat=as.matrix(priors[,3:ncol(priors)]) ; rownames(priorMat)=priors[,1] ; colnames(priorMat)=colnames(priors)[3:ncol(priors)]
out=setup() ##5 objects
priorWeights=out[[1]]#PriorWeights: rows are genes, cols are TFs
gene.train=out[[4]];gene.test=out[[5]]##gene.train and gene.test rows are obs and cols are genes
tf.train=out[[2]];tf.test=out[[3]]##TF train and TF test: rows are obs and cols are TFs
geneNames=as.character(gene.exp[,2]) ##gene names


##BART in random noise
sigma_noise=1
y=rnorm(nrow(tf.train),mean=0,sd=sigma_noise)
##no thin
bart.mod=bart(x.train=tf.train,y.train=y,ntree=10,nskip=3000,ndpost=5000)
sum(bart.mod$varcount)
prop_calc(bart.mod)
max(prop_calc(bart.mod))
which.max(prop_calc(bart.mod))

##Thin
bart.mod=bart(x.train=tf.train,y.train=y,ntree=10,nskip=3000,ndpost=5000,keepevery=100)
sum(bart.mod$varcount)
prop_calc(bart.mod)
max(prop_calc(bart.mod))
which.max(prop_calc(bart.mod))


out=list()
nsim=50
count=1
sig_list=c(.1,.5,1,2,4,8)
for(sigma in sig_list){
  out[[count]]=list()
  out[[count]][["sigma"]]=sigma
  maxvec=numeric(nsim)
  maxvec_thin=numeric(nsim)
  maxvec_superthin=numeric(nsim)
  for(i in 1:nsim){
    y=rnorm(nrow(tf.train),mean=0,sd=sigma)
    ##no thin
    bart.mod=bart(x.train=tf.train,y.train=y,ntree=10,nskip=2000,ndpost=5000,verbose=F)
    maxvec[i]=max(prop_calc(bart.mod))
    ##thin
    bart.mod=bart(x.train=tf.train,y.train=y,ntree=10,nskip=2000,ndpost=5000,keepevery=25,verbose=F)
    maxvec_thin[i]=max(prop_calc(bart.mod))
    ##super thin
    bart.mod=bart(x.train=tf.train,y.train=y,ntree=10,nskip=2000,ndpost=5000,keepevery=100,verbose=F)
    maxvec_superthin[i]=max(prop_calc(bart.mod))
    if(i%%10==0) print(i)
  }
  out[[count]][["maxvec"]]=maxvec
  out[[count]][["maxvec_thin"]]=maxvec_thin
  out[[count]][["maxvec_super"]]=maxvec_superthin
  count=count+1
}

summary(out[[1]][["maxvec"]])
summary(out[[1]][["maxvec_thin"]])
summary(out[[1]][["maxvec_thin"]])


