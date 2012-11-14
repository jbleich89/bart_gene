##bart multicore main

# setwd("~/Research_Genomics/bart_gene") ##work
# setwd("~/Documents/Research/Genomics/bart_gene/") ##home
source("bart_mc_fns.R")

library(BayesTree)
library(multicore)
library(parallel)
##work
#priors=read.table("~/Research_Genomics/CHIP.priorprobs.39.txt",header=T) ##work
#gene.exp=read.table("~/Research_Genomics/expression.genes.txt",header=T) ##work
#tf.exp=read.table("~/Research_Genomics/expression.tfs.39.txt",header=T) ##work
##Home
#priors=read.table("~/Documents/Research/Genomics/Real_Data/CHIP.priorprobs.39.txt",header=T) ##home
#gene.exp=read.table("~/Documents/Research/Genomics/Real_Data/expression.genes.txt",header=T) ##home
#tf.exp=read.table("~/Documents/Research/Genomics/Real_Data/expression.tfs.39.txt",header=T) ##home
##variance
priors=read.table("~/bart_genome/CHIP.priorprobs.39.txt",header=T) ##MC
gene.exp=read.table("~/bart_genome/expression.genes.txt",header=T) ##MC
tf.exp=read.table("~/bart_genome/expression.tfs.39.txt",header=T) ##MC

num_cores=22
##objects-priorWeights (cols), priorMat (probs), gene.train, tf.train, geneNames

priorMat=as.matrix(priors[,3:ncol(priors)]) ; rownames(priorMat)=priors[,1] ; colnames(priorMat)=colnames(priors)[3:ncol(priors)]
out=setup() ##5 objects
priorWeights=out[[1]]#PriorWeights: rows are genes, cols are TFs
gene.train=out[[4]];gene.test=out[[5]]##gene.train and gene.test rows are obs and cols are genes
tf.train=out[[2]];tf.test=out[[3]]##TF train and TF test: rows are obs and cols are TFs
geneNames=as.character(gene.exp[,2]) ##gene names


run_BART=function(geneList,tf.mat,runBoot=F,nboot=100,num_cores=(detectCores()-1),...){
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
    #bart.mod=bart(tf.mat.prior,gene.response,ntree=10,nskip=2000,ndpost=5000,keepevery=25)
    bart.mod=bart(x.train=tf.mat.prior,y.train=gene.response,,...) ## run BART
    strung.reps=rep(names(reps),times=as.integer(reps)) ##get names of TFs by number of reps
    var_prop=prop_calc_prior(bart.mod,strung.reps) ##calculate var counts using prior
    sumVars=sum(sum_calc_prior(bart.mod,strung.reps)) ##total sum of variables
    out[[gene]][["var"]]=var_prop ##var props
    out[[gene]][["yhat"]]=(bart.mod$yhat.train.mean) ##in sample yhat
    out[[gene]][["sig"]]=mean(bart.mod$sigma) ##estimate of sigma
    out[[gene]][["numSplits"]]=sumVars ##total number of splits used in whole model
    count=count+1
    print("bart done")
    if(runBoot==T) {#runs non-par boot
      
      # part 1- return boot matrix
      ##use multiple cores here
      print("getting boot")
      boot_list=mclapply(1:nboot,getBootIter,gene.vec=gene.response,tf.train=tf.mat,mc.cores=num_cores,...)
      boot_mat=t(do.call(cbind,boot_list))
      colnames(boot_mat)=colnames(tf.mat)
      print("boot done")
      #out[[gene]][["numSplitsNull"]]=boot_mat_list[["nullSum"]]
      #boot_mat=getBootMat(gene.vec=gene.response,tf.train=tf.train,strung.rep.vec=strung.reps,ntree=5)
      ##part 2 do computations
      if(count<=10) out[[gene]][["boot_mat"]]=boot_mat
      
      ##Simult. Coverage-Bands
      boot.se=apply(boot_mat,2,sd)
      mean_boot=apply(boot_mat,2,mean)
      coverConst=bisectK(tol=.01,coverage=.95,boot_mat=boot_mat,x_left=1,x_right=20,countLimit=100)
      out[[gene]][["const"]]=coverConst
      print(coverConst)
      mean(sapply(1:nrow(boot_mat), function(s) all(boot_mat[s,]-mean_boot<=coverConst*boot.se)))
      simul_trueTFs=var_prop[which(var_prop>=mean_boot+coverConst*boot.se)]
      out[[gene]][["s_trueTF"]]=simul_trueTFs
      ##Simult. Coverage-Max
      maxcut=quantile(apply(boot_mat,1,max),.95)
      maxtrue_TFs=var_prop[which(var_prop>=maxcut)]
      out[[gene]][["max_trueTF"]]=max_trueTFs
      ##pointwise coverage
      q95_point=apply(boot_mat,2,quantile, probs=.95)
      point_trueTFs=var_prop[which(var_prop>=q95_point)]
      out[[gene]][["p_trueTF"]]=point_trueTFs
      ##leave everything else in this call
      print("choosing TFs done")
    }
    
    if(count%%10==0) print(count)
  }
  t1=Sys.time()
  print(t1-t0)
  return(out)
}

##run barts
#out5=run_BART(geneList=geneNames[1:100],tf.train,runBoot=T,ntree=5,nskip=1000,ndpost=2000,verbose=F)
#save(out5,file="out5_1.rdata")
out10=run_BART(geneList=geneNames[1:5],tf.train,runBoot=T,num_cores=num_cores,ntree=10,nskip=1000,ndpost=2000,verbose=F)
#save(out5,file="out10_1.rdata")
#out20=run_BART(geneList=geneNames[1:100],tf.train,runBoot=F,ntree=20,nskip=1000,ndpost=2000,verbose=F)
#save(out5,file="out20_1.rdata")







