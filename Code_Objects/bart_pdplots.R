##Partial Dependence Plots Work
require(BayesTree)
setwd("~/Documents/Research/Genomics")
source("bart_fns.R")
source("bart_fns2.R")

##Set-up
set.seed(20)

##TF settings
tf.size1=c(10,15,20,25)
tf.size2=c(10,20,30,40) ##script is set to work with this right now


##Observation settings
n50=50
n100=100
n300=300

mean.tf=0 #1.4862e-05 ##generate the X matrix which will be fixed- sl
sd.tf=1 #.4778817-not using just to keep larger numbers at play

##Generate Design Matrices

tf.exp.300=sapply(rep(n300,max(tf.size2)),rnorm,mean=mean.tf,sd=sd.tf) #gives full matrix
tf.exp.50=tf.exp.300[1:n50,]
tf.exp.100=tf.exp.300[1:n100,]

##Beta settings

tf.beta.1=c(1,rep(0,times=max(tf.size2-1)))
#tf.beta.2=c(2,rep(0,times=max(tf.size2-1))) not really using here


##Function params
n.tree.vec=c(5,10,15,20)
factor.vec1=c(0,.25,.5,1,1.5,2,3.5,5,8)
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12) ##work with this right now 
burn=2000
post=15000


get_pd_plot=function(tf.exp,n.obs,tf.beta,n.tree,factor,tf.size,burn_size,post_size,
                     sigest=NA,thin.size=1,plot_vars){
  signal.temp=sum(abs(tf.exp%*%tf.beta))/n.obs
  sigma=signal.temp*factor
  gene.exp=as.numeric(tf.exp%*%tf.beta+rnorm(n.obs,mean=0,sd=sigma))
  train.exp=tf.exp[1:n.obs,1:tf.size]
  ##For case of N-use sample SD-else if will use whatever is pre-specified
  if(n.obs<=ncol(train.exp) & is.na(sigest)) bart.sig=sd(gene.exp)
  else bart.sig=NA
  print(dim(train.exp))
  print(bart.sig)
  bart.mod = pdbart(x.train=train.exp,
                  y.train=gene.exp,
                  xind=plot_vars,
                  ntree=n.tree.vec,
                  sigest=bart.sig,
                  nskip=burn_size,
                  ndpost=post_size,
                  keepevery=thin.size,
                  verbose=T) ##Suppress printing
  return(bart.mod)
}


##3 Different Settings of Factor Levels
mod2=get_pd_plot(tf.exp=tf.exp.300,n.obs=300,tf.beta=tf.beta.1,n.tree=20,
            factor=1,tf.size=40,burn_size=burn,post_size=post,thin.size=5,plot_vars=c(1,2))

mod2=get_pd_plot(tf.exp=tf.exp.300,n.obs=300,tf.beta=tf.beta.1,n.tree=20,
                 factor=5,tf.size=40,burn_size=burn,post_size=post,thin.size=5,plot_vars=c(1,2))

mod2=get_pd_plot(tf.exp=tf.exp.300,n.obs=300,tf.beta=tf.beta.1,n.tree=20,
                 factor=7,tf.size=40,burn_size=burn,post_size=post,thin.size=5,plot_vars=c(1,2))
names(mod2)
mod2$xlbs
par(mfrow=c(1,2))
plot(mod2,main="Partial Dependence Plot\n 5x N/S")
