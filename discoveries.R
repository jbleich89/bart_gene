load("out10t_500g.rdata")
length(out10)
genes=geneNames[1:500]
apply(out10[[genes[3]]]$boot_mat,2,quantile,probs=.95)



out10[[genes[3]]][["max_trueTF"]]
out10[[genes[3]]][["p_trueTF"]]
out10[[genes[3]]][["s_trueTF"]]
out10[[genes[4]]][["var"]]
quantile(apply(out10[[genes[3]]][["boot_mat"]],1,max),.95)

maxcounts=numeric(ncol(tf.train))
names(maxcounts)=colnames(tf.train)
name=names(maxcounts)[15]
name
names(maxcounts)
genes=geneNames[3]
y=getGeneResponse(genes)
t=getGeneTestResponse(genes)
fb=bart(x.train=tf.train,x.test=tf.test,y.train=y,ntree=200,nskip=500,ndpost=1000)
rb=bart(x.train=(tf.train[,names(out10[[genes[3]]][["s_trueTF"]])]),x.test=(tf.test[,names(out10[[genes[3]]][["s_trueTF"]])]),y.train=y,ntree=200,nskip=500,ndpost=1000)
1-sum((t-fb$yhat.test.mean)^2)/sum(t^2)
1-sum((t-rb$yhat.test.mean)^2)/sum(t^2)
sum(t^2)
sum((t)^2)
sum((t-mean(t))^2)
length(t)
summary(lm(y~))



(getDiscoveryCounts(counts=maxcounts,genes=geneNames[1:500],"max_trueTF"))
getDiscoveryProps(counts=maxcounts,genes=geneNames[1:500],"s_trueTF")
getDiscoveryProps(counts=maxcounts,genes=geneNames[1:500],"max_trueTF")

a=getDiscoveryCounts(counts=maxcounts,genes=geneNames[1:500],"p_trueTF")
b=getDiscoveryCounts(counts=maxcounts,genes=geneNames[1:500],"s_trueTF")
c=getDiscoveryCounts(counts=maxcounts,genes=geneNames[1:500],"max_trueTF")

##code for bar plots 
comb=rbind(a,b-a)

barplot(a,las=3,,ylab="Inclusion Count",main="Inclusion Counts for TFs")
barplot(b,las=3,,ylab="Inclusion Count",main="Inclusion Counts for TFs",xlim)
barplot(c,las=3,,ylab="Inclusion Count",main="Inclusion Counts for TFs")
barplot(comb,beside=F,las=3,ylab="Inclusion Count",main="Inclusion Counts for TFs",col=c("blue","red")
        ,cex.names=.90,axis.lty=1)
##



##functions
getDiscoveryProps=function(counts,genes,idx){
  for(name in names(counts)){
    tf_temp=0
    for(gene in genes){
      temp=name %in% names(out10[[gene]][[idx]])
      tf_temp=tf_temp+temp
    }
    counts[which(names(counts)==name)]=tf_temp
  }  
  return(counts/length(genes))
}

getDiscoveryCounts=function(counts,genes,idx){
  for(name in names(counts)){
    tf_temp=0
    for(gene in genes){
      temp=name %in% names(out10[[gene]][[idx]])
      tf_temp=tf_temp+temp
    }
    counts[which(names(counts)==name)]=tf_temp
  }  
  return(counts)
}

maxcounts
name

source("C:/Users/jbleich/Desktop/Dropbox/BART_gene/load_genes.R")
transcription_to_gene
sum(do.call(cbind,lapply(transcription_to_gene,length)))