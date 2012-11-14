##bart multicore functions
#Functions for real data


##Functions###########
##Rounding function for setting prior
roundUpNice <- function(x, nice=c(1,2,3,4,5,6,7,8,9,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

changePrior=function(location){
  temp=location*10
  ifelse(temp>1,roundUpNice(temp),1)
}

ChangePriorVec=function(locVec){
  sapply(locVec,changePrior)
}
####################
##Functions for calculating proportions
sum_calc=function(bart.obj){apply(bart.obj$varcount,2,sum)}

##calculate variable proportions 
prop_calc=function(bart.obj){
  tot=sum(apply(bart.obj$varcount,2,sum))
  props=apply(bart.obj$varcount,2,sum)/tot
  round(props,3)
}

#1.
##calculate variable sum
sum_calc_prior=function(bart.obj,prior.vec){
  temp=apply(bart.obj$varcount,2,sum)
  return(tapply(temp,prior.vec,sum))
}

#2.
##calculate variable proportions 
prop_calc_prior=function(bart.obj,prior.vec){
  temp=sum_calc_prior(bart.obj,prior.vec)
  props=temp/sum(temp)
  return(round(props,3))
}

####################################
#Bootstrap functions



################
#BART Loop Functions

getGeneResponse=function(gene){
  idx=which(colnames(gene.train)==gene)
  gv=as.numeric(gene.train[,idx])
  names(gv)=rownames(gene.train)
  return(gv)
}

getReps=function(gene){
  return(priorWeights[gene,])
}



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


##bootstrap for multicore
getbootIter=function(iter,gene.vec,tf.train,...)
  iter ##hack just to make sure mc works
  n=length(gene.vec)
  perm.sample=gene.vec[sample(1:n,n,F)] ##permute y vector
  tf.mat.null=tf.train
  bart.boot=bart(x.train=tf.mat.null,y.train=perm.sample,...)
  props=prop_calc(bart.obj=bart.boot)
  return(props )
##regular bootstrap function


############################
##Setup function
##This function creates:
#PriorWeights: rows are genes, cols are TFs
##gene.train and gene.test rows are obs and cols are genes
##TF train and TF test: rows are obs and cols are TFs
setup=function(){
  ##remove first two columns
  prior_nums=priors[,3:ncol(priors)]
  geneNames=as.character(gene.exp[,2])
  
  ##get weights- rownames are genes. colnames are TFs
  priorWeights=apply(prior_nums,2,ChangePriorVec)
  rownames(priorWeights)=geneNames
  
  ##prior col 1 goes with gene col 2 and vice versa
  priorWeights[,1]==gene.exp[,2]
  
  
  ##need transpose of tf mat
  tf.names=tf.exp[,1]
  temp.tf=tf.exp[,5:ncol(tf.exp)]
  tf.full=t(temp.tf) ##This is training matrix
  colnames(tf.full)=tf.names
  rownames(tf.full)=colnames(tf.exp)[5:ncol(tf.exp)]
  dim(tf.full) #check 314 conditions and 39 TFs
  hold=ceiling(.2*nrow(tf.full))
  
  #get training and test
  train.end=314-hold
  tf.train=tf.full[1:train.end,]
  tf.test=tf.full[(train.end+1):nrow(tf.full),]
  
  ##Clean gene matrix
  geneMatClean=t(gene.exp[,5:ncol(gene.exp)])
  rownames(geneMatClean)=rownames(tf.full)
  colnames(geneMatClean)=geneNames
  
  ##Get gene test and train
  gene.train=geneMatClean[1:train.end,]
  gene.test=geneMatClean[(1+train.end):nrow(geneMatClean),]
  return(list(priorWeights,tf.train,tf.test,gene.train,gene.test))
}