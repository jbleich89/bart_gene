setwd("~/Documents/Research/Genomics/Real_Data")

priors=read.table("CHIP.priorprobs.39.txt",header=T)
gene.exp=read.table("expression.genes.txt",header=T)
tf.exp=read.table("expression.tfs.39.txt",header=T)
dim(gene.exp)
dim(tf.exp)
head(gene.exp)

head(priors)

##remove first two columns
prior_nums=priors[,3:ncol(priors)]
##get weights
weightMat=apply(prior_nums,2,ChangePriorVec)
#reattach gene info
priorWeights=data.frame(priors[,1:2],weightMat)
##prior col 1 goes with gene col 2 and vice versa
priorWeights[,1]==gene.exp[,2]
head(gene.exp)

##need transpose of tf mat
tf.names=tf.exp[,1]
temp.tf=tf.exp[,5:ncol(tf.exp)]
tf.train=t(temp.tf) ##This is training matrix
colnames(tf.train)=tf.names
##do names line up
colnames(tf.train)==colnames(priorWeights[3:41])

sort(table(gene.exp[,1]),descending=T)
which(gene.exp[,1]=="HXT12")
gene.exp[2846:2847,1:5]
##Functions###########3
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

pigamma()
paste