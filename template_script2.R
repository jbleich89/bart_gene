setwd("~/Research_Genomics/bart_gene") ##work
source("real_functions.R")
library(BayesTree)


priors=read.table("~/Research_Genomics/CHIP.priorprobs.39.txt",header=T) ##work
gene.exp=read.table("~/Research_Genomics/expression.genes.txt",header=T) ##work
tf.exp=read.table("~/Research_Genomics/expression.tfs.39.txt",header=T) ##work


out=setup() ##5 objects
priorWeights=out[[1]]#PriorWeights: rows are genes, cols are TFs
gene.train=out[[4]];gene.test=out[[5]]##gene.train and gene.test rows are obs and cols are genes
tf.train=out[[2]];tf.test=out[[3]]##TF train and TF test: rows are obs and cols are TFs
geneNames=as.character(gene.exp[,2]) ##gene names

run_BART(geneList=gene,tf.train,runBoot=T,ntree=10,nskip=1000,ndpost=2000)

run_BART=function(geneList,tf.mat,runBoot=F,...){
  out=list()
  count=0
  for(gene in geneList){
    out[[gene]]=list()
    out[[gene]][["name"]]=gene #store name
    
    ##Setup train
    reps=getReps(gene) ##get reps
    tf.mat.prior=tf.train[,rep(1:ncol(tf.mat),reps)] ##make longer matrix
    print(dim(tf.mat.prior))
    
    gene.response=getGeneResponse(gene) ##get response
    length(gene.response) ##should be 251
    #bart.mod=bart(tf.mat.prior,gene.response,ntree=10,nskip=2500,ndpost=2500)
    bart.mod=bart(x.train=tf.mat.prior,y.train=gene.response,,...) ## run BART
    strung.reps=rep(names(reps),times=as.integer(reps)) ##get names of TFs by number of reps
    var_prop=prop_calc_prior(bart.mod,strung.reps) ##calculate var counts using prior
    out[[gene]][["var"]]=var_prop ##var props
    out[[gene]][["yhat"]]=(bart.mod$yhat.train.mean) ##in sample yhat
    out[[gene]][["sig"]]=mean(bart.mod$sigma) ##estimate of sigma
    count=count+1
    if(runBoot==T) nonpar_boot(gene.vec=gene.response,strung.rep.vec=strung.reps,tf.train=tf.mat,count=count,,...) #runs non-par boot
    ##Lets unwrap the function for the most part
    # part 1- return boot matrix
    ##part 2 do computations
    ##leave everything else in this call
    
    
    
    if(count%%10==0) print(count)
  }
  return(out)
}

gene=geneNames[3]
tf.mat=tf.train
plot(bart.mod)
var_prop
nboot=100
gene.vec=gene.response
strung.rep.vec=strung.reps


nonpar_boot=function(gene.vec,strung.rep.vec,tf.train,nboot=10,count,...){ ##this procedure adjusts for simult.
  ##needs original training matrix
  print(count)
  n=length(gene.vec) ##for permutation
  boot_mat=matrix(0,nrow=nboot,ncol=ncol(priorWeights)) ##ncol is fixed
  for(j in 1:nboot){ 
    perm.sample=gene.vec[sample(1:n,n,F)] ##permute y vector
    print(perm.sample[1])
    tf.mat.break=priorBreak(tf.train,length(strung.rep.vec))
    #bart.boot=bart(x.train=tf.mat.break,y.train=perm.sample,ntree=20,nskip=2000,ndpost=2000,verbose=F)
    bart.boot=bart(x.train=tf.mat.break,y.train=perm.sample,...)
    props=prop_calc_prior(bart.boot,strung.rep.vec)
    boot_mat[j,]=props ##need max to adjust for simult.
    if(j%%25==0) print(j)
  }
  if(count<=10) out[[gene]][["boot_mat"]]=boot_mat 
  boot.se=apply(boot_mat,2,sd)
  mean_boot=apply(boot_mat,2,mean)
  coverConst=bisectK(tol=.01,coverage=.95,boot_mat=boot_mat,x_left=1,x_right=20,limit=100)
  mean(sapply(1:nboot, function(s) all(boot_mat[s,]-mean_boot<=coverConst*boot.se)))
  trueTFs=which(var_prop>=mean_boot+coverConst*boot.se)
  return(trueTFs)
}


##function for breaking prior
priorBreak=function(tf.mat,lengthPrior,numTFs=39){
  extraTFs=lengthPrior-numTFs
  samp=sample(1:numTFs,extraTFs,replace=T) ##uniform sample with replacement from other TFs 
  extraCols=tf.mat[,samp]
  return(cbind(tf.mat,extraCols))
}


##bisection method for finding simult. coverage 
bisectK=function(tol,coverage,boot_mat,x_left,x_right,limit){
  count=0
  x_left=2
  x_right=15
  guess=(x_left+x_right)/2
  while(.5*(x_right-x_left)>=tol & count<countLimit){
    mean_boot=apply(boot_mat,2,mean)
    empCoverage=mean(sapply(1:nboot, function(s) all(boot_mat[s,]-mean_boot<=guess*boot.se)))
    if(empCoverage-coverage==0) break
    else if((empCoverage-coverage)<0) x_left=guess
    else x_right=guess
    guess=(x_left+x_right)/2
    count=count+1
  }
  return(guess)
}






#####################OTHER STUFF FOR NOW####################

##Run some models
bart5=run_BART(geneNameVec[1:1000],ntree=5,ndpost=3000,nskip=2000,verbose=F)
bart10=run_BART(geneNameVec[1:1000],ntree=10,ndpost=3000,nskip=2000,verbose=F)
bart20=run_BART(geneNameVec[1:1000],ntree=20,ndpost=3000,nskip=2000,verbose=F)


cor5=sapply(1:1000, function(x) cor(bart5[[x]][[2]],priorWeights[x,]))
cor10=sapply(1:1000, function(x) cor(bart10[[x]][[2]],priorWeights[x,]))
cor20=sapply(1:1000, function(x) cor(bart20[[x]][[2]],priorWeights[x,]))

summary(cor5,na.rm=T)
summary(cor10,na.rm=T)
summary(cor20,na.rm=T)

max5=sapply(1:1000, getMaxes, bartList=bart5)
max10=sapply(1:1000, getMaxes, bartList=bart10)
max20=sapply(1:1000, getMaxes, bartList=bart20)
which(max20==max(max20))
summary(max5)
summary(max10)
summary(max20)
Sbart10[[732]][[1]]
priorWeights["YDL056W",]

##
c=max(bart5[[1]][[2]])
d=which(bart5[[1]][[2]]==c)
names(d)
idx=1
bartList=bart5
getMaxes=function(bartList,idx){
  mx=max(bartList[[idx]][[2]])
  nameMx=which(bartList[[idx]][[2]]==mx)
  names(mx)=names(nameMx)[1]
  return(mx)
}
names()
cor(x[[2]][[2]],priorWeights[2,])
which(x[[3]][[2]]==max(x[[3]][[2]]))
x[[3]][[2]]
reps
  ##Need to get prior function stuff from computer 
prioy=getReps(gene)
priorWeights[gene,]
train.exp=train.exp.temp[,rep(1:ncol(train.exp.temp),reps)]




head(gene.exp)
gene.exp[3,5:318]==gene.response



