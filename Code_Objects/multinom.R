##Sampling max of multinomial
numTFs=39
pTF=1/numTFs
probVec=rep(pTF,times=numTFs)
N=24

rmultinom

res=rmultinom(10000,N,prob=probVec)
res[,1]
sum(bart.boot$varcount[,ncol(bart.boot$varcount)])
1+2+2+2+3+2+2+3+3+4
maxVec=apply(res,2,max)
quantile(maxVec,.95)
mean(apply(bart.boot$varcount,1,sum))

avg_prop=function(bart.mod){
vc=bart.boot$varcount
sum(apply(vc,2,sum))
t=apply(vc,2,mean)  
sum(t)
round(t/sum(t),3)
mean(apply(vc,2,sum))
vc[1,]/sum(vc[1,])
  
}
prop_calc(bart.boot)


dim(vc)

res=numeric(nrow(vc))
for(i in 1:50){
  x=rmultinom(50,16,probVec)
  x[1,]
  r=apply(x,1,sum)/sum(apply(x,1,sum))
  max(r)
  res[i]=
}

