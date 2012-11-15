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


nrows=nrow(tf.train)
ncols=ncol(tf.train)
x.train=matrix(rnorm(nrows*ncols,0,1),nrow=nrows,ncol=ncols)

out=list()
nsim=100
count=1
sig_list=c(.1,.5,1,2,4,8)
start=Sys.time()
for(sigma in sig_list){
  out[[count]]=list()
  out[[count]][["sigma"]]=sigma
  maxvec=numeric(nsim)
  maxvec_thin=numeric(nsim)
  maxvec_superthin=numeric(nsim)
  for(i in 1:nsim){
    y=rnorm(nrow(tf.train),mean=0,sd=sigma)
    ##no thin
    bart.mod=bart(x.train=x.train,y.train=y,ntree=10,nskip=2000,ndpost=5000,verbose=F)
    maxvec[i]=max(prop_calc(bart.mod))
    ##thin
    bart.mod=bart(x.train=x.train,y.train=y,ntree=10,nskip=2000,ndpost=5000,keepevery=25,verbose=F)
    maxvec_thin[i]=max(prop_calc(bart.mod))
    ##super thin
    bart.mod=bart(x.train=x.train,y.train=y,ntree=10,nskip=2000,ndpost=5000,keepevery=100,verbose=F)
    maxvec_superthin[i]=max(prop_calc(bart.mod))
    if(i%%10==0) print(i)
  }
  out[[count]][["maxvec"]]=maxvec
  out[[count]][["maxvec_thin"]]=maxvec_thin
  out[[count]][["maxvec_super"]]=maxvec_superthin
  count=count+1
}
print(Sys.time()-start)
hist(out[[3]][["maxvec"]],col="grey",main="Max Incl. Freqs. Across Datasets",xlab="Max Inclusion Freq.",breaks=30)
hist(out[[3]][["maxvec_thin"]],col="grey",main="Max Incl. Freqs. Across Datasets-Thinned",xlab="Max Inclusion Freq.",breaks=30)

hist(out[[6]][["maxvec"]],col="grey",main="Max Incl. Freqs. Across Datasets-",xlab="Max Inclusion Freq.",breaks=30)
hist(out[[6]][["maxvec_thin"]],col="grey",main="Max Incl. Freqs. Across Datasets-Thinned",xlab="Max Inclusion Freq.",breaks=30)



summary(out[[1]][["maxvec_super"]])




##
out=list()
par(mgp=c(1.8,.5,0),mar=c(3,3,2,1))
par(mfrow=c(1,2))
nsim=50
ymin=0

nrows=nrow(tf.train)
ncols=ncol(tf.train)
x.train=matrix(rnorm(nrows*ncols,0,1),nrow=nrows,ncol=ncols)
sig_list=c(.1,.5,1,2,4,8)
sigma=sig_list[3]
sigma
n=12000

y=rnorm(nrows,0,sigma)
bart.noise=bart(x.train=x.train,y.train=y,ntree=10,nskip=1,ndpost=n,keepevery=1)
prop_calc(bart.noise)
cum_props[n,]

dim(cum_props)
cums=apply(bart.noise$varcount,2,cumsum)
cum_props=cums/apply(cums,1,sum)
ymax=max(cum_props[2000:n,])
ymax
rain=rainbow(ncol(tf.train))

plot(1:n,cum_props[,1],type="l",col=rain[1],ylim=c(ymin,ymax),main="Cumulative Inclusion Props.",
     xlab="Gibbs Sample",ylab="Inclusion Freq.")
sapply(2:ncol(tf.train), function(s) points(1:n,cum_props[,s],type="l",col=rain[s]))
##follow max
maxs=apply(cum_props,1,max)
mins=apply(cum_props,1,min)
plot(1:n,maxs,ylim=c(ymin,ymax),type="l",col="red",lwd=3, main="Max and Min Cumulative Props.",ylab="Inclusion Freq.",xlab="Gibbs Sample")
points(1:n,mins,type="l",col="blue",lwd=3)



##multiple chains
b1=bart(x.train=tf.train,y.train=y,ntree=10,nskip=1,ndpost=n/4,keepevery=1)
b2=bart(x.train=tf.train,y.train=y,ntree=10,nskip=1,ndpost=n/4,keepevery=1)
b3=bart(x.train=tf.train,y.train=y,ntree=10,nskip=1,ndpost=n/4,keepevery=1)
b4=bart(x.train=tf.train,y.train=y,ntree=10,nskip=1,ndpost=n/4,keepevery=1)
fullvar=rbind(b1$varcount,b2$varcount,b3$varcount,b4$varcount)
max(prop_calc(b4))
which.max(prop_calc(b4))

fullcums=apply(fullvar,2,cumsum)
fullcum_props=fullcums/apply(fullcums,1,sum)
rain=rainbow(ncol(tf.train))
dim(fullcums)

plot(1:n,fullcum_props[,1],type="l",col=rain[1],ylim=c(min(fullcum_props),max(fullcum_props[1000:n,])),main="Cumulative Inclusion Props.",
     xlab="Gibbs Sample",ylab="Inclusion Freq.")
sapply(2:ncol(tf.train), function(s) points(1:n,fullcum_props[,s],type="l",col=rain[s]))
##follow max
maxs=apply(fullcum_props,1,max)
mins=apply(fullcum_props,1,min)
plot(1:n,maxs,ylim=c(0,.2),type="l",col="red",lwd=3, main="Max and Min Cumulative Props.",ylab="Inclusion Freq.",xlab="Gibbs Sample")
points(1:n,mins,type="l",col="blue",lwd=3)

