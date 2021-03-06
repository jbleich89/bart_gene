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

##put in setup
priorMat=as.matrix(priors[,3:ncol(priors)])
rownames(priorMat)=priors[,1] ; colnames(priorMat)=colnames(priors)[3:ncol(priors)]
out=setup() ##5 objects
priorWeights=out[[1]]#PriorWeights: rows are genes, cols are TFs
gene.train=out[[4]];gene.test=out[[5]]##gene.train and gene.test rows are obs and cols are genes
tf.train=out[[2]];tf.test=out[[3]]##TF train and TF test: rows are obs and cols are TFs
geneNames=as.character(gene.exp[,2]) ##gene names


##run barts
out5=run_BART(geneList=geneNames[1:100],tf.train,runBoot=T,ntree=5,nskip=1000,ndpost=2000,verbose=F)
#save(out5,file="out5_1.rdata")
out10=run_BART(geneList=geneNames[1:100],tf.train,runBoot=F,ntree=10,nskip=1000,ndpost=2000,verbose=F)
#save(out5,file="out10_1.rdata")
out20=run_BART(geneList=geneNames[1:100],tf.train,runBoot=F,ntree=20,nskip=1000,ndpost=2000,verbose=F)
#save(out5,file="out20_1.rdata")


out20=run_BART(geneList=geneNames[3],tf.train,runBoot=F,ntree=20,nskip=1000,ndpost=2000,verbose=F)


genelist=geneNames[3]
genelist
tf.mat=tf.train
gene=geneNames[10]
plot(1:ncol(tf.mat),var_prop,pch=16)


##Consider FDR?
run_BART=function(geneList,tf.mat,runBoot=F,...){
  t0=Sys.time() ##timer
  out=list() ##set up list
  count=0
  for(gene in geneList){ ##for each gene you're considering
    out[[gene]]=list() ##make a list
    out[[gene]][["name"]]=gene #store name
    
    ##Setup train
    reps=getReps(gene) ##get reps. gets the right prior weights 
    tf.mat.prior=tf.train[,rep(1:ncol(tf.mat),reps)] ##make longer matrix with row reps
    print(dim(tf.mat.prior)) ##check for right dimension
    
    gene.response=getGeneResponse(gene) ##get response- returns it from big matrix of response
    length(gene.response) ##should be 251
    bart.mod=bart(tf.mat.prior,gene.response,ntree=10,nskip=2000,ndpost=5000,keepevery=25)
    #bart.mod=bart(x.train=tf.mat.prior,y.train=gene.response,,...) ## run BART
    strung.reps=rep(names(reps),times=as.integer(reps)) ##get names of TFs by number of reps
    var_prop=prop_calc_prior(bart.mod,strung.reps) ##calculate var counts using prior
    sumVars=sum(sum_calc_prior(bart.mod,strung.reps)) ##total sum of variables
    out[[gene]][["var"]]=var_prop ##var props
    out[[gene]][["yhat"]]=(bart.mod$yhat.train.mean) ##in sample yhat
    out[[gene]][["sig"]]=mean(bart.mod$sigma) ##estimate of sigma
    out[[gene]][["numSplits"]]=sumVars ##total number of splits used in whole model
    count=count+1
    if(runBoot==T) {#runs non-par boot
      
      # part 1- return boot matrix
      boot_mat_list=getBootMat(gene.vec=gene.response,tf.train=tf.train,strung.rep.vec=strung.reps,...) ##generate boot matrix
      boot_mat=boot_mat_list[["boot"]] ##gets
      out[[gene]][["numSplitsNull"]]=boot_mat_list[["nullSum"]]
      #boot_mat=getBootMat(gene.vec=gene.response,tf.train=tf.train,strung.rep.vec=strung.reps,ntree=5)
      ##part 2 do computations
      if(count<=10) out[[gene]][["boot_mat"]]=boot_mat
  
      ##Simult. Coverage Part
      boot.se=apply(boot_mat,2,sd)
      mean_boot=apply(boot_mat,2,mean)
      coverConst=bisectK(tol=.01,coverage=.95,boot_mat=boot_mat,x_left=1,x_right=20,countLimit=100)
      out[[gene]][["const"]]=coverConst
      print(coverConst)
      mean(sapply(1:nrow(boot_mat), function(s) all(boot_mat[s,]-mean_boot<=coverConst*boot.se)))
      simul_trueTFs=which(var_prop>=mean_boot+coverConst*boot.se)
      out[[gene]][["s_trueTF"]]=simul_trueTFs
      
      ##pointwise coverage
      q95_point=apply(boot_mat,2,quantile, probs=.95)
      point_trueTFs=which(var_prop>=q95_point)
      out[[gene]][["p_trueTF"]]=point_trueTFs
      ##leave everything else in this call
    }
    
      if(count%%10==0) print(count)
  }
  t1=Sys.time()
  print(t1-t0)
  return(out)
}

var_prop
maxid=apply(boot_mat,1,which.max)
maxs=apply(boot_mat,1,max)
maxcut=quantile(maxs,.95)
table(maxid)
hist(maxs)
maxcut
mean_boot=apply(boot_mat,2,mean); boot.se=apply(boot_mat,2,sd)
coverConst=bisectK(tol=.1,coverage=.95,boot_mat=boot_mat,x_left=1,x_right=20,countLimit=100)
coverConst

mean(sapply(1:nrow(boot_mat), function(s) all(boot_mat[s,]-mean_boot<=coverConst*boot.se)))
simul_trueTFs=which(var_prop>=mean_boot+coverConst*boot.se)
simul_trueTFs
plot(1:ncol(tf.mat),var_prop,pch=16,xlab="TF",main=paste("Variable Inclusion Plot for",gene,"\n Simultaneous"))
abline(h=maxcut,col="red")
sapply(1:ncol(tf.mat),function(s) segments(s,0,x1=s,mean_boot[s]+coverConst*boot.se[s],col="blue"))
sapply(1:ncol(tf.mat),function(s) segments(s,0,x1=s,quantile(boot_mat[,s],.95),col="red"))
par(mgp=c(1.8,.5,0))

###
getBootMat=function(gene.vec,tf.train,strung.rep.vec,...,nboot=10){ ##needs original training matrix for simult. inference
  out=list()
  n=length(gene.vec) ##for permutation
  boot_mat=matrix(0,nrow=nboot,ncol=ncol(tf.train)) ##ncol is fixed-at 39  
  nullSum=numeric(nboot) ##vector of total vars used
  for(j in 1:nboot){ 
    perm.sample=gene.vec[sample(1:n,n,F)] ##permute y vector
    #print(perm.sample[1])
    tf.mat.null=tf.train  ##use training matrix with no priors
    bart.boot=bart(x.train=tf.mat.null,y.train=perm.sample,ntree=10,nskip=2000,ndpost=5000,keepevery=25)
    #bart.boot=bart(x.train=tf.mat.null,y.train=perm.sample,...)
    props=prop_calc(bart.obj=bart.boot) ##CHECK! THIS COULD BE WRONG. Think its ok 11/4
    nullSum[j]=sum(sum_calc(bart.obj=bart.boot)) ##total number of splits. 
    #print(length(props))
    boot_mat[j,]=props ##need max to adjust for simult.
    if(j%%25==0) print(j)
  }
  out[["boot"]]=boot_mat
  out[["nullSum"]]=nullSum
  return(out)
}


##function for breaking prior
##Now Defunct
# priorBreak=function(tf.mat,lengthPrior,numTFs=39){
#   extraTFs=lengthPrior-numTFs
#   samp=sample(1:numTFs,extraTFs,replace=T) ##uniform sample with replacement from other TFs 
#   extraCols=tf.mat[,samp]
#   temp=(cbind(tf.mat,extraCols))
#   colnames(temp)=c(colnames(tf.mat),colnames(tf.mat)[samp])
#   return(temp)
# }


##bisection method for finding simult. coverage 
##now defunct
bisectK=function(tol,coverage,boot_mat,x_left,x_right,countLimit){
  count=0
  x_left=x_left
  x_right=x_right
  guess=(x_left+x_right)/2
  while(.5*(x_right-x_left)>=tol & count<countLimit){
    mean_boot=apply(boot_mat,2,mean)
    boot.se=apply(boot_mat,2,sd)
    empCoverage=mean(sapply(1:nrow(boot_mat), function(s) all(boot_mat[s,]-mean_boot<=guess*boot.se)))
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

train.exp=train.exp.temp[,rep(1:ncol(train.exp.temp),reps)]




head(gene.exp)
gene.exp[3,5:318]==gene.response



