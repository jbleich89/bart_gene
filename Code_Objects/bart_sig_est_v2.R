require(BayesTree)
setwd("~/Documents/Research/Genomics")
source("bart_fns.R")
source("bart_fns2.R")

##Set-up
set.seed(20)

##TF settings
tf.size2=c(10,20,30,40) ##script is set to work with this right now

##Observation settings
n50=50
n300=300

mean.tf=0 #1.4862e-05 ##generate the X matrix which will be fixed- sl
sd.tf=1 #.4778817-not using just to keep larger numbers at play

##Generate Design Matrices

tf.exp.300=sapply(rep(n300,max(tf.size2)),rnorm,mean=mean.tf,sd=sd.tf) #gives full matrix
tf.exp.50=tf.exp.300[1:n50,]
tf.exp.100=tf.exp.300[1:n100,]

##Beta settings

tf.beta.1=c(1,rep(0,times=max(tf.size2-1)))
#tf.beta.2=c(2,rep(0,times=max(tf.size2-1))) not really using here


##Function params
n.tree.vec=c(10,20)
factor.vec2=c(1,2,3,4,5) ##work with this right now 
burn=2000
post=5000
n.genes=100

##generate list of genes
gene.list2=get_gene_data_list(n.obs=50,tf.exp.50,factor.vec2,tf.beta.1,reps=n.genes)

##Run Algos
bart.sig.50.prior=bart_sig_prior(tf.exp.50,gene.list2,n.obs=50,reps=n.genes,factor.vec2,n.tree.vec,tf.vec=40,
             sigest=NA,burn_size=burn,post_size=post,
                                 sig_df=3,sig_quant=.9)

ls.sig.est=numeric(length(bart.sig.50.prior))
for(i in 1:length(bart.sig.50.prior)){
  ls.sig.est[i]=bart.sig.50.prior[[i]][[3]]
}

prior.sig.vec=numeric(length(bart.sig.50.prior))
for(i in 1:length(bart.sig.50.prior)){
  prior.sig.vec[i]=bart.sig.50.prior[[i]][[2]][4]
}
prior.sig.vec
ls.sig.est

##how many times in total number of reps was right TF selected
prop.true.prior=numeric(length(bart.sig.50.prior))
for(i in 1:length(bart.sig.50.prior)){
  prop.true.prior[i]=bart.sig.50.prior[[i]][[4]]
}
prop.true.prior

##Avg. variable usage proportion
prop.selected.prior=numeric(length(bart.sig.50.prior))
for(i in 1:length(bart.sig.50.prior)){
  prop.selected.prior[i]=bart.sig.50.prior[[i]][[5]]
}
prop.selected.prior


##use large df and center at median
bart.sig.50.post=bart_sig_post(tf.exp.50,gene.list2,n50,reps=n.genes,factor.vec2,n.tree.vec,tf.vec=40,
                               sigest=prior.sig.vec,burn_size=burn,post_size=post,
                               sig_df=200,sig_quant=.5)

post.sig.vec=numeric(length(bart.sig.50.post))
for(i in 1:length(bart.sig.50.post)){
  post.sig.vec[i]=bart.sig.50.post[[i]][[2]][4]
}
post.sig.vec


prop.true.post=numeric(length(bart.sig.50.post))
for(i in 1:length(bart.sig.50.post)){
  prop.true.post[i]=bart.sig.50.post[[i]][[4]]
}
prop.true.post



prop.selected.post=numeric(length(bart.sig.50.post))
for(i in 1:length(bart.sig.50.post)){
  prop.selected.post[i]=bart.sig.50.post[[i]][[5]]
}
prop.selected.post

##Default settings for df and quantile, just using estimated sigma
bart.sig.50.post.def=bart_sig_post(tf.exp.50,gene.list2,n50,reps=n.genes,factor.vec2,n.tree.vec,tf.vec=40,
                               sigest=prior.sig.vec,burn_size=burn,post_size=post,
                               sig_df=3,sig_quant=.9)

post.sig.vec.def=numeric(length(bart.sig.50.post.def))
for(i in 1:length(bart.sig.50.post.def)){
  post.sig.vec.def[i]=bart.sig.50.post.def[[i]][[2]][4]
}
post.sig.vec.def


prop.true.post.def=numeric(length(bart.sig.50.post.def))
for(i in 1:length(bart.sig.50.post.def)){
  prop.true.post.def[i]=bart.sig.50.post.def[[i]][[4]]
}
prop.true.post.def


prop.selected.post.def=numeric(length(bart.sig.50.post.def))
for(i in 1:length(bart.sig.50.post.def)){
  prop.selected.post.def[i]=bart.sig.50.post.def[[i]][[5]]
}
prop.selected.post.def

##Bind together to see desired outcomes

##Generate column and row names
cnames=c("Prior_Prop. ","Post_Prop._Fixed ", "Post_Prop._Default")
##rownames
rnames=character(length(n.tree.vec)*length(factor.vec2))
count=1
for(i in 1:length(n.tree.vec)){
  for(j in 1:length(factor.vec2)){
    rnames[count]=paste(n.tree.vec[i]," Trees ", factor.vec2[j], "x N/S", sep="")
    count=count+1
  }
}

##Proportion of n=rep simulations where true TF was selected
propmat=cbind(prop.true.prior,prop.true.post,prop.true.post.def)
colnames(propmat)=cnames
rownames(propmat)=rnames
propmat

##Avg. selection proportion for True TF across simulations
selmat=cbind(prop.selected.prior,prop.selected.post,prop.selected.post.def)
colnames(selmat)=cnames
rownames(selmat)=rnames
round(selmat,3)


##Compare 3 variances
sig.mat=round(cbind(ls.sig.est,prior.sig.vec,post.sig.vec,post.sig.vec.def),2)
colnames(sig.mat)=c("LS_Est. ","Prior_Est. ","Post_Est._Fixed ","Post_Est._Default")
rownames(sig.mat)=rnames
sig.mat


par(mfrow=c(1,2))
plot_sig_full(bart.sig.50.post.def[[2]][[1]],n_burn=1500,subtitle="Post_Default")
plot_sig_full(bart.sig.50.post[[2]][[1]],n_burn=1500,subtitle="Posterior_Fixed")


mean(bart.sig.50.prior[[1]][[1]]$sigma)
bart.sig.50.prior[[2]][[1]]$sigest
mean(bart.sig.50.prior[[2]][[1]]$sigma)
bart.sig.50.post[[2]][[1]]$sigma[1:5]
plot(bart.sig.50.prior[[2]][[1]])
bart.sig.50.prior[[2]][[1]]$first.sigma
bart.sig.50.prior[[1]][[2]][1]


bart_sig_summary(bart.sig.50.prior)
bart_sig_summary(bart.sig.50.post)



save(bart.sig.50.prior,file="prior50.R")
save(bart.sig.50.post,file="post50.R")

load("prior50.R")
load("post50.R")

##Functions

##1.
##Creates y-vectors- first index in list corresponds to noise factor right now. 
##Can change to additional element to serve as legend or hashmap it. 
get_gene_data_list=function(n.obs,tf.exp,factor.vec, tf.beta,reps){
  signal.temp=sum(abs(tf.exp%*%tf.beta))/n.obs
  out=list()
  for(i in 1:length(factor.vec)){
    out[[i]]=list()
    sigma=signal.temp*factor.vec[i]
    for(j in 1:reps){
      out[[i]][[j]]=as.numeric(tf.exp%*%tf.beta+rnorm(n.obs,mean=0,sd=sigma))
    }
  }
  return(out)
}


##2.
bart_sig_prior=function(tf.exp,gene_list,n.obs,reps,factor.vec,tree.vec,
                        tf.vec,sigest=NA,burn_size,post_size,sig_df=3,sig_quant=.9,thin.size=2){
##Second item in the list has the following order:
  #(tree size, number of tfs, noise factor, posterior sigma estimate)
  big.list=list()
  bart.sig=sigest
  count=1
  for(i in 1:length(tree.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(tf.vec)){
        sig.holder=numeric(reps)
        prior.sig=numeric(reps)
        true.tf=numeric(reps)
        sel.prop=numeric(reps) ##vector for proportion of time selected
        ##general loop-creates sigma, then response, then runs BART
        for(rep in 1:reps){
          gene.exp=gene_list[[j]][[rep]]
          #print(gene.exp[1:4])
          print(var(gene.exp))
          train.exp=tf.exp[1:n.obs,1:tf.vec[k]]
  
          bart.mod = bart(x.train=train.exp,
                          y.train=gene.exp,
                          ntree=n.tree.vec[k],
                          sigest=bart.sig,
                          sigdf=sig_df,
                          sigquant=sig_quant,
                          nskip=burn_size,
                          ndpost=post_size,
                          keepevery=thin.size,
                          verbose=F) ##Suppress printing
          print(rep)
          sig.holder[rep]=mean(bart.mod$sigma) ##hold the mean of the posterior-estimate of sigma
          prior.sig[rep]=bart.mod$sigest ##hold the priors for comparison 
          ##add component to check how many times it was correct
          sums=sum_calc(bart.mod)
          true.tf[rep]=which.max(sums)==1
          sel.prop[rep]=sums[1]/sum(sums)
        }
        ##printing
        print(mean(sig.holder))
        print(bart.mod$sigest)
        print(true.tf)
        
        big.list[[count]]=list()
        big.list[[count]][[1]]=bart.mod
        big.list[[count]][[2]]=c(tree.vec[i],tf.vec[k],factor.vec[j],mean(sig.holder))
        big.list[[count]][[3]]=mean(prior.sig)
        big.list[[count]][[4]]=mean(true.tf)
        big.list[[count]][[5]]=mean(sel.prop)
        print(count)
        count=count+1
      }
    }
  }
  return(big.list)
}


##3.
bart_sig_post=function(tf.exp,gene_list,n.obs,reps,factor.vec,tree.vec,
                       tf.vec,sigest.vec,burn_size,post_size,sig_df=3,sig_quant=.90,thin.size=2){
  ##Second item in the list has the following order:
  #(tree size, number of tfs, noise factor, posterior sigma estimate)
  ##Need to make more robust to handle different orderings-right now just copies above
  big.list=list()
  count=1
  for(i in 1:length(tree.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(tf.vec)){
        sig.holder=numeric(reps) ##posterior var vector
        true.tf=numeric(reps) ##vector for correct proportion
        sel.prop=numeric(reps) ##vector for proportion of time selected
        for(rep in 1:reps){
          
          gene.exp=gene_list[[j]][[rep]]
          #print(gene.exp[1:4])
          print(var(gene.exp))
          
          train.exp=tf.exp[1:n.obs,1:tf.vec[k]]
          bart.mod = bart(x.train=train.exp,
                          y.train=gene.exp,
                          ntree=n.tree.vec[k],
                          sigest=sigest.vec[count],
                          sigdf=sig_df,
                          sigquant=sig_quant,
                          nskip=burn_size,
                          ndpost=post_size,
                          keepevery=thin.size,
                          verbose=F) ##Suppress printing
          print(rep)
          sig.holder[rep]=mean(bart.mod$sigma)
          sums=sum_calc(bart.mod)
          true.tf[rep]=which.max(sums)==1
          sel.prop[rep]=sums[1]/sum(sums)
        }
        ##printing to check
        print(mean(sig.holder))
        print(bart.mod$sigest)
        print(true.tf)
        
        big.list[[count]]=list()
        big.list[[count]][[1]]=bart.mod
        big.list[[count]][[2]]=c(tree.vec[i],tf.vec[k],factor.vec[j],mean(sig.holder))
        big.list[[count]][[3]]=NA
        big.list[[count]][[4]]=mean(true.tf)
        big.list[[count]][[5]]=mean(sel.prop)
        
        print(count)
        count=count+1
      }
    }
  }
  return(big.list)
}

##4. 
plot_sig_full=function(bart.list,n_burn,subtitle){
  sig=c(bart.list$first.sigma,bart.list$sigma)
  burn.idx=1:n_burn
  burn=sig[burn.idx]
  post=sig[-burn.idx]
  plot(burn,main=paste("Convergence Plot for Sigma\n",subtitle),
       xlab="Sample Number",
       ylab="Sigma",
       xlim=c(0,length(sig)),
       ylim=c(min(sig),max(sig)),
       col="red"
  )
  points((length(burn.idx)+1):length(sig),post)
}

##5. 
bart_sig_summary=function(bart.mod){
  output=data.frame()
  count=1
  for(i in 1:length(bart.mod)){
    output[count,1]=bart.mod[[i]][[2]][1] ##Trees
    output[count,2]=bart.mod[[i]][[2]][2] ##TF
    output[count,3]=bart.mod[[i]][[2]][3] ##Noise
    
    sums=sum_calc(bart.mod[[i]][[1]])
    props=prop_calc(bart.mod[[i]][[1]])
    sorted.sums=sort(sums,decreasing=T)
    sorted.props=sort(props,decreasing=T)    
    
    output[count,4]=which.max(sums)==1  ##Does our appear most?               
    output[count,5]=sorted.props[1] ##Top three TFs
    output[count,6]=which(props==sorted.props[1])[1] ##How to deal with ties?
    output[count,7]=sorted.props[2]
    output[count,8]=which(props==sorted.props[2])[1]
    output[count,9]=sorted.props[3]
    output[count,10]=which(props==sorted.props[3])[1]
    output[count,11]=sums[1]
    output[count,12]=sum(sums) ##Total splits 
    output[count,13]=mean(bart.mod[[i]][[1]]$sigma)
    output[count,14]=props[1]
    print(count)
    count=count+1
  }
  colnames(output)=c("Num_Trees","Num_TFs","Noise/Signal","True_Most_Common?","1st","Name-1st",
                     "2nd","Name-2nd","3rd","Name-3rd","True_TF_Sum","Tot_Splits","Sigma_Estimate","True_TF_Prop")
  return(output)
}


##simulations with chisquare 
library(geoR)
d=150
x=seq(.5,3,.01)
dens=dinvchisq(x,df=d,scale=1)
plot(x,dens)
n=rinvchisq(10000,d,scale=1)
mean(n)
quantile(n,c(.25,.5,.75))
plot(n)
var(n)