setwd("~/Research_Genomics/bart_gene") ##work
source("real_functions.R")



priors=read.table("~/Research_Genomics/CHIP.priorprobs.39.txt",header=T) ##work
gene.exp=read.table("~/Research_Genomics/expression.genes.txt",header=T) ##work
tf.exp=read.table("~/Research_Genomics/expression.tfs.39.txt",header=T) ##work



##remove first two columns
prior_nums=priors[,3:ncol(priors)]

##get weights
weightMat=apply(prior_nums,2,ChangePriorVec)

#reattach gene info
priorWeights=data.frame(priors[,1:2],weightMat)

##prior col 1 goes with gene col 2 and vice versa
priorWeights[,1]==gene.exp[,2]


##need transpose of tf mat
tf.names=tf.exp[,1]
temp.tf=tf.exp[,5:ncol(tf.exp)]
tf.train=t(temp.tf) ##This is training matrix
colnames(tf.train)=tf.names
##do names line up
colnames(tf.train)==colnames(priorWeights[3:41])

geneNameVec=gene.exp[,2] ##version has 6026 genes

prepForBart=function(gene)

for(gene in geneNameVec) 
  
  
  ##Need to get prior function stuff from computer 











