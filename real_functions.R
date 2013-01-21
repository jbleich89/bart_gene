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
  
#  ##get weights- rownames are genes. colnames are TFs
#  priorWeights=apply(prior_nums,2,ChangePriorVec)
#  rownames(priorWeights)=geneNames
#  
#  ##prior col 1 goes with gene col 2 and vice versa
#  priorWeights[,1]==gene.exp[,2]
  
  
  ##need transpose of tf mat
  tf.names=tf.exp[,1]
  temp.tf=tf.exp[,5:ncol(tf.exp)]
  tf.full=t(temp.tf) ##This is training matrix
  colnames(tf.full)=tf.names
  rownames(tf.full)=colnames(tf.exp)[5:ncol(tf.exp)]
  dim(tf.full) #check 314 conditions and 39 TFs
  
  ##randomize later
  cv_hold=ceiling(.1*nrow(tf.full))
  test_hold=ceiling(.1*nrow(tf.full))
		  
  #get training and test
  train.end=nrow(tf.full)-(cv_hold + test_hold)
  tf.train=tf.full[1:train.end,]
  tf.cv=tf.full[(train.end+1):(train.end+cv_hold),]
  tf.test=tf.full[(train.end+cv_hold+1):nrow(tf.full),]
  
  ##Clean gene matrix
  geneMatClean=t(gene.exp[,5:ncol(gene.exp)])
  rownames(geneMatClean)=rownames(tf.full)
  colnames(geneMatClean)=geneNames
  
  ##Get gene test and train
  gene.train=geneMatClean[1:train.end,]
  gene.cv=geneMatClean[(train.end+1):(train.end+cv_hold),]
  gene.test=geneMatClean[(train.end+cv_hold+1):nrow(geneMatClean),]
  list(
	  tf.train=tf.train,
	  tf.cv=tf.cv, 
	  tf.test=tf.test,
	  gene.train=gene.train, 
	  gene.cv=gene.cv, 
	  gene.test=gene.test
  )
}