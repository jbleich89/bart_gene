require(BayesTree)
setwd("~/Documents/Research/Genomics")
setwd("~/Research_Genomics/bart_gene/Code_Objects")
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
num.tf=40
tf.beta.1=c(1,rep(0,times=max(tf.size)-1))
tf.beta.3=c(1,.5,.5,rep(0,times=num.tf-3))
#tf.beta.2=c(2,rep(0,times=max(tf.size2-1))) not really using here


##Function params
n.tree.vec=c(5,10,15,20)
factor.vec1=c(0,.25,.5,1,1.5,2,3.5,5,8)
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12) ##work with this right now 
burn=2500
post=5000

##Prior weight vectors
rep2=c(rep(1,times=2),2:max(tf.size2))
rep5=c(rep(1,times=5),2:max(tf.size2))
rep10=c(rep(1,times=10),2:max(tf.size2))
rep50=c(rep(1,times=50),2:max(tf.size2))

##for tf 3
rep121=c(1,2,2,3:num.tf)

prior.vec=rep121
tf.size=num.tf
##Set up number of repititions and repeats in data frame
rep.nums=sapply(1:tf.size,function(x) length(which(prior.vec==x)))
#print(rep.nums)
reps=rep.nums

#print(reps)
dim(train.exp)
tf.exp=tf.exp.300
tf.beta=tf.beta.3
n=nrow(tf.exp)
sigma=2*sum(abs(tf.exp%*%tf.beta))/n
train.exp=tf.exp[,rep(1:ncol(tf.exp),reps)]
gene=as.numeric(tf.exp%*%tf.beta+rnorm(n,mean=0,sd=sigma))

tf.null=tf.exp

bart.true=bart(x.train=train.exp,y.train=gene,ntree=10,nskip=2000,ndpost=5000,keepevery=25)
var_prop=prop_calc_prior(bart.true,rep121)

boot=100
boot_mat=matrix(0,nrow=boot,ncol=ncol(tf.null))
sums=numeric(boot)
for(i in 1:boot){
  perm=gene[sample(1:n,n,replace=F)]
  bart.false=bart(x.train=tf.null,y.train=perm,ntree=10,nskip=1000,ndpost=5000,keepevery=25,verbose=F)
  boot_mat[i,]=prop_calc(bart.false)
  sums[i]=sum(bart.false$varcount)
     print(i)
}

boot_mat[1,]
maxid=apply(boot_mat,1,which.max)
maxs=apply(boot_mat,1,max)
maxcut=quantile(maxs,.95)
table(maxid)
hist(maxs)

summary(sums)
mean_boot=apply(boot_mat,2,mean); boot.se=apply(boot_mat,2,sd)
coverConst=bisectK(tol=.1,coverage=.95,boot_mat=boot_mat,x_left=1,x_right=20,countLimit=100)
coverConst
mean(sapply(1:nrow(boot_mat), function(s) all(boot_mat[s,]-mean_boot<=coverConst*boot.se)))
simul_trueTFs=which(var_prop>=mean_boot+coverConst*boot.se)
simul_trueTFs
plot(1:40,var_prop,main="Simulated TFs at Noise of 2x\n Simultaneous",xlab="TF")
points(1:3,var_prop[1:3],pch=16)
abline(h=maxcut,col="red")

sapply(1:ncol(tf.exp),function(s) segments(s,0,x1=s,mean_boot[s]+coverConst*boot.se[s],col="blue"))
sapply(1:ncol(tf.exp),function(s) segments(s,0,x1=s,quantile(boot_mat[,s],.95),col="red"))