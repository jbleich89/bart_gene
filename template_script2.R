setwd("~/Research_Genomics/bart_gene") ##work
source("real_functions.R")
library(BayesTree)


priors=read.table("~/Research_Genomics/CHIP.priorprobs.39.txt",header=T) ##work
gene.exp=read.table("~/Research_Genomics/expression.genes.txt",header=T) ##work
tf.exp=read.table("~/Research_Genomics/expression.tfs.39.txt",header=T) ##work



##remove first two columns
prior_nums=priors[,3:ncol(priors)]
geneNames=as.character(gene.exp[,2])

##get weights
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

##do names line up
colnames(tf.train)==colnames(priorWeights[3:41])

##Set up response
head(gene.exp)


geneNameVec=gene.exp[,2] ##version has 6026 genes

prepForBart=function(gene){
  idx=which(priorWeights[,2]==gene)
  w.vec=(priorWeights[gene,])  
}


run_BART=function(geneList,...){
  out=list()
  count=0
  for(gene in geneList){
    out[[gene]]=list()
    out[[gene]][[1]]=gene #store name
    
    ##Setup train
    reps=getReps(gene) ##get reps
    tf.train.prior=tf.train[,rep(1:ncol(tf.train),reps)]
    print(dim(tf.train.prior))
    
    gene.response=getGeneResponse(gene)
    length(gene.response)
    bart.mod=bart(x.train=tf.train.prior,y.train=gene.response,...) ## run BART
    strung.reps=rep(names(reps),times=as.integer(reps))
    var_prop=prop_calc_prior(bart.mod,strung.reps)
    out[[gene]][[2]]=var_prop
    out[[gene]][[3]]=(bart.mod$yhat.train.mean)
    out[[gene]][[4]]=mean(bart.mod$sigma)
    count=count+1
    if(count%%10==0) print(count)
  }
  return(out)
}

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

getGeneResponse=function(gene){
  idx=which(colnames(gene.train)==gene)
  gv=as.numeric(gene.train[,idx])
  names(gv)=rownames(gene.train)
  return(gv)
}

getReps=function(gene){
  return(priorWeights[gene,])
}


head(gene.exp)
gene.exp[3,5:318]==gene.response



