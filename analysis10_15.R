##Initial Analysis:

#save(out10,file="tree10_100.rdata")
#save(out5,file="tree5_100.rdata")

#load("tree5_100.rdata")
#load("tree10_100.rdata")

gene=geneNames[1]
ngenes=100

##Count number of true TFs
## pointwise
p_count5=0 ; p_count10=0 ; s_count5=0 ; s_count10=0

for(i in 1:ngenes){
p_count5= p_count5+length(out5[[geneNames[i]]][["p_trueTF"]])
p_count10= p_count10+length(out10[[geneNames[i]]][["p_trueTF"]])
s_count5= s_count5+length(out5[[geneNames[i]]][["s_trueTF"]])
s_count10= s_count10+length(out10[[geneNames[i]]][["s_trueTF"]])
}
p_count5
p_count10
s_count5
s_count10

##Correlations for Ntrees
corrList=numeric(ngenes)
for(i in 1:ngenes){
corrList[i]=cor(out5[[geneNames[i]]][["var"]],out10[[geneNames[i]]][["var"]])
}
summary(corrList)
hist(corrList,main="Correlations of Incl.Freqs \n for m=5 and m=10",col="grey",breaks=25)

cor(priorMat,priorWeights)
class(priorMat)
((x-mean(x)/sd(x))%*%((y-mean(y))/sd(y)))/250


##Corrs of TF exp and Var incl freq
cor5=sapply(1:ngenes, function(x) cor(out5[[geneNames[x]]][["var"]],priorMat[x,]))
cor10=sapply(1:ngenes, function(x) cor(out10[[geneNames[x]]][["var"]],priorMat[x,]))
hist(cor5,col="grey",breaks=30,main="Correlations for 5 Trees")
hist(cor10,col="grey",breaks=35,main="Correlations for 10 Trees")


na.rows= which(is.na(cor10)) ##same for 10


##Shane's Plots

#5 Trees
lik.cor5=numeric(ngenes)
for(i in 1:ngenes){
  cor.temp=apply(tf.train,2, function(x) cor(gene.train[,i],x))
  lik.cor5[i]=cor(cor.temp,out5[[geneNames[i]]][["var"]])
}
lowlim=min(min(lik.cor5),min(cor5,na.rm=T))
uplim=max(max(lik.cor5),max(cor5,na.rm=T))
plot(cor5,lik.cor5, xlim=c(-.3,.7), ylim=c(-.3,.7),pch=16,xlab="Prior and Incl. Freq. Corr.",
     ylab="Corr. of Data Corrs. and Incl. Freq.", main="Likelihood vs. Prior Correlations \n 5 Trees")
abline(0,1,col="red",lwd=3)


##10 Trees
lik.cor10=numeric(ngenes)
for(i in 1:ngenes){
  cor.temp=apply(tf.train,2, function(x) cor(gene.train[,i],x))
  lik.cor10[i]=cor(cor.temp,out10[[geneNames[i]]][["var"]])
}


lowlim=min(min(lik.cor10),min(cor10,na.rm=T))
uplim=max(max(lik.cor10),max(cor10,na.rm=T))
plot(cor10,lik.cor10, xlim=c(lowlim,uplim), ylim=c(lowlim,uplim),pch=16,xlab="Prior and Incl. Freq. Corr.",
     ylab="Corr. of Data Corrs. and Incl. Freq.", main="Likelihood vs. Prior Correlations \n 10 Trees")
abline(0,1,col="red",lwd=3)


  