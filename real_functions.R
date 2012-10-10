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

bart.obj=bart.mod
prior.vec=reps
length(temp)
length(reps)
t=rep(reps,reps)
t=rep(names(reps),times=as.integer(reps))
names(reps)
length(t)

nonpar_boot=function(gene.vec,strung.rep.vec,nboot=100){
  n=length(gene.vec)
  maxVec=numeric(nboot)
  for(j in 1:nboot){ ##this procedure adjusts for simult.
  perm.sample=gene.vec[sample(1:n,n,F)]
  bart.boot=bart(x.train=tf.train.prior,y.train=perm.sample,...)
  props=prop_calc_prior(bart.boot,strung.rep.vec)
  maxVec[j]=max(props) ##need max to adjust for simult.
  }
  return(maxVec)
}