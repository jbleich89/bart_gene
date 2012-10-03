require(BayesTree)
setwd("~/Documents/Research/Genomics")
source("bart_fns.R")
source("bart_fns2.R")

##Set-up
set.seed(20)

##TF settings
tf.size2=c(10,20,30,40) ##script is set to work with this right now

##Observation settings
n15=15
n50=50
n300=300

mean.tf=0 #1.4862e-05 ##generate the X matrix which will be fixed- sl
sd.tf=1 #.4778817-not using just to keep larger numbers at play

##Generate Design Matrices

tf.exp.300=sapply(rep(n300,max(tf.size2)),rnorm,mean=mean.tf,sd=sd.tf) #gives full matrix
tf.exp.50=tf.exp.300[1:n50,]
tf.exp.100=tf.exp.300[1:n100,]
tf.exp.15=tf.exp.300[1:n15,]
##Beta settings

tf.beta.1=c(1,rep(0,times=max(tf.size2-1)))
#tf.beta.2=c(2,rep(0,times=max(tf.size2-1))) not really using here


##Function params
n.tree.vec=c(10,20)
factor.vec.sig=c(1,2,3,4,5) ##work with this right now 
burn=2000
post=4000
n.genes=30
n.obs=10
tf.exp.input=tf.exp.300[1:n.obs,]
signal.input=sum(abs(tf.exp.input%*%tf.beta.1))/nrow(tf.exp.input)
truth=tf.exp.input%*%tf.beta.1

##generate list of genes
gene.list2=get_gene_data_list(n.obs=n.obs,tf.exp.input,factor.vec.sig,tf.beta.1,reps=n.genes)

test.list.sig=gen_test(n.obs=20,num.cols=40, tf.beta=tf.beta.1,factor.vec=factor.vec.sig,signal=signal.input)

##Run Algos
bart.sig.50.prior=bart_sig_prior_test(tf.exp.input,gene.list2,n.obs=n.obs,reps=n.genes,factor.vec.sig,n.tree.vec,tf.vec=40,
             sigest=NA,burn_size=burn,post_size=post,
                                 sig_df=3,sig_quant=.9,test.list=test.list,truth=truth)

ls.sig.est=sapply(1:length(bart.sig.50.prior),function(i) bart.sig.50.prior[[i]][[3]])
prior.sig.vec=sapply(1:length(bart.sig.50.prior),function(i) bart.sig.50.prior[[i]][[2]][4])
##how many times in total number of reps was right TF selected
prop.true.prior=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[4]])
##Avg. variable usage proportion
prop.selected.prior=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[5]])

####
##use large df and center at median
bart.sig.50.post=bart_sig_post_test(tf.exp.input,gene.list2,n.obs=n.obs,reps=n.genes,factor.vec.sig,n.tree.vec,tf.vec=40,
                               sigest.vec=prior.sig.vec,burn_size=burn,post_size=post,
                               sig_df=200,sig_quant=.5,test.list=test.list, truth=truth)

post.sig.vec=sapply(1:length(bart.sig.50.post),function(i) bart.sig.50.post[[i]][[2]][4])
##how many times in total number of reps was right TF selected
prop.true.post=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post[[i]][[4]])
##Avg. variable usage proportion
prop.selected.post=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post[[i]][[5]])


##Default settings for df and quantile, just using estimated sigma
bart.sig.50.post.def=bart_sig_post_test(tf.exp.input,gene.list2,n.obs=n.obs,reps=n.genes,factor.vec.sig,n.tree.vec,tf.vec=40,
                               sigest.vec=prior.sig.vec,burn_size=burn,post_size=post,
                               sig_df=3,sig_quant=.9,test.list=test.list, truth=truth)

post.sig.vec.def=sapply(1:length(bart.sig.50.post.def),function(i) bart.sig.50.post.def[[i]][[2]][4])
##how many times in total number of reps was right TF selected
prop.true.post.def=sapply(1:length(bart.sig.50.post.def), function(i) bart.sig.50.post.def[[i]][[4]])
##Avg. variable usage proportion
prop.selected.post.def=sapply(1:length(bart.sig.50.post.def),
                              function(i) bart.sig.50.post.def[[i]][[5]])

##Bind together to see desired outcomes

##Generate column and row names
cnames=c("Prior_Prop. ","Post_Prop._Fixed ", "Post_Prop._Default")
##rownames
rnames=character(length(n.tree.vec)*length(factor.vec.sig))
count=1
for(i in 1:length(n.tree.vec)){
  for(j in 1:length(factor.vec.sig)){
    rnames[count]=paste(n.tree.vec[i]," Trees ", factor.vec.sig[j], "x N/S", sep="")
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
sigma.vec=rep(signal.input*(factor.vec.sig),times=2)
sig.mat=round(cbind(sigma.vec,ls.sig.est,prior.sig.vec,post.sig.vec,post.sig.vec.def),2)
colnames(sig.mat)=c("Sigma","LS_Est. ","Prior_Est. ","Post_Est._Fixed ","Post_Est._Default")
rownames(sig.mat)=rnames
sig.mat

##Predictive Performance:
##prior
#prior.perf.train=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[6]] )
#prior.perf.test=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[7]] )
prior.truth.rmse=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[8]] )
##centered post
#post.perf.train=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post[[i]][[6]] )
#post.perf.test=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post[[i]][[7]] )
post.truth.rmse=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post[[i]][[8]] )
##def. post
#post.d.perf.train=sapply(1:length(bart.sig.50.post.def), function(i) bart.sig.50.post.def[[i]][[6]] )
#post.d.perf.test=sapply(1:length(bart.sig.50.post.def), function(i) bart.sig.50.post.def[[i]][[7]] )
post.d.truth.rmse=sapply(1:length(bart.sig.50.post.def), function(i) bart.sig.50.post.def[[i]][[8]] )

##Training matrix
#train.rmse.mat=round(cbind(prior.perf.train,post.perf.train,post.d.perf.train),3)
#test.rmse.mat=round(cbind(prior.perf.test,post.perf.test,post.d.perf.test),3)
truth.rmse.mat=(cbind(prior.truth.rmse,post.truth.rmse,post.d.truth.rmse))
#train.rmse.mat
#test.rmse.mat
truth.rmse.mat


##Coverage
prior.cov=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[9]] )
post.cov=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post.def[[i]][[9]] )
post.d.cov=sapply(1:length(bart.sig.50.post.def), function(i) bart.sig.50.post.def[[i]][[9]] )
cbind(prior.cov,post.cov,post.d.cov)

#Width 90%
prior.wid=sapply(1:length(bart.sig.50.prior), function(i) bart.sig.50.prior[[i]][[10]] )
post.wid=sapply(1:length(bart.sig.50.post), function(i) bart.sig.50.post.def[[i]][[10]] )
post.d.wid=sapply(1:length(bart.sig.50.post.def), function(i) bart.sig.50.post.def[[i]][[10]] )
cbind(prior.wid,post.wid,post.d.wid)

x=apply(bart.sig.50.prior[[2]][[1]]$yhat.train, 2, quantile, probs=c(.05,.95))
y=apply(bart.sig.50.post[[2]][[1]]$yhat.train, 2, quantile, probs=c(.05,.95))
t=tf.exp.input%*%tf.beta.1
sum(sapply(1:length(truth), function(i) truth[i]>=y[1,i] & truth[i]<=y[2,i]))/length(t)

mean(sapply(1:ncol(y),function(i) abs(y[1,i]-y[2,i])))


temp=(bart.sig.50.prior[[2]][[1]]$sigma)
yrange=c(min(temp),quantile(temp,.99))
par(mfrow=c(1,2))
plot_sig_full(bart.sig.50.prior[[2]][[1]],n_burn=1000,subtitle="Initial Run",y.lims=yrange)
plot_sig_full(bart.sig.50.post[[2]][[1]],n_burn=1000,subtitle="Shared Variance",yrange)

bart.sig.50.prior[[1]][[8]]

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
get_gene_data_list=function(n.obs,tf.exp,factor.vec,tf.beta,reps){
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

##Get test data
gen_test=function(n.obs=250,num.cols=40,tf.beta,factor.vec,signal){
  out=list()
  mean.tf=0 #1.4862e-05 ##generate the X matrix which will be fixed- sl
  sd.tf=1 #.4778817-not using just to keep larger numbers at play
  test.mat=sapply(rep(n.obs,num.cols),rnorm,mean=mean.tf,sd=sd.tf)
  for(i in 1:length(factor.vec)){
    
    out[[i]]=list()
    sigma=signal*factor.vec[i]
    out[[i]][[1]]=as.numeric(test.mat%*%tf.beta+rnorm(n.obs,mean=0,sd=sigma))
    out[[i]][[2]]=factor.vec[i]
  }
  out[[length(factor.vec)+1]]=list()
  out[[length(factor.vec)+1]][[1]]=test.mat 
  out[[length(factor.vec)+1]][[2]]=NA
  return(out)
}

##2.
bart_sig_prior_test=function(tf.exp,gene_list,n.obs,reps,factor.vec,tree.vec,
                        tf.vec,sigest=NA,burn_size,post_size,sig_df=3,
                        sig_quant=.9,thin.size=2,test.list,truth){
##Second item in the list has the following order:
  #(tree size, number of tfs, noise factor, posterior sigma estimate)
  big.list=list()
  bart.sig=sigest
  count=1
  test_factor_idx=sapply(1:length(test.list), function(i) test.list[[i]][[2]])
  test.mat=test.list[[which(is.na(test_factor_idx))]][[1]] ##get test mat
  
  for(i in 1:length(tree.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(tf.vec)){
        
        sig.holder=numeric(reps)
        prior.sig=numeric(reps)
        true.tf=numeric(reps)
        sel.prop=numeric(reps) ##vector for proportion of time selected
        
        train.rmse=numeric(reps) ##in sample RMSE
        test.rmse=numeric(reps) ##vector for RMSE on test set
        true.rmse=numeric(reps) ##vector for RMSE on truth
        coverage=numeric(reps)  
        width=numeric(reps)
        ##general loop-creates sigma, then response, then runs BART
        for(rep in 1:reps){
          gene.exp=gene_list[[j]][[rep]]
          #print(gene.exp[1:4])
          
          train.exp=tf.exp[1:n.obs,1:tf.vec[k]]
          
          ##Handle n less than p          
          if(n.obs<=ncol(train.exp) & is.na(sigest)) bart.sig=sd(gene.exp)
          else bart.sig=NA
          print(bart.sig)
  
          bart.mod = bart(x.train=train.exp,
                          y.train=gene.exp,
                          #x.test=test.mat,
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
          ##test info
         # train.rmse[reps]=sqrt(sum((bart.mod$yhat.train.mean-gene.exp)^2)/n.obs)
          #loc=which(test_factor_idx==factor.vec[j])
          #y.test=test.list[[loc]][[1]]
          #test.rmse[reps]=sqrt(sum((bart.mod$yhat.test.mean-y.test)^2)/nrow(test.mat))
          #train.rmse[reps]=sqrt(sum((bart.mod$yhat.train.mean-gene.exp)^2)/n.obs)
          true.rmse[rep]=sqrt(sum((bart.mod$yhat.train.mean-truth)^2)/nrow(train.exp))
          
          ##coverage
          cov.temp=apply(bart.mod$yhat.train,2,quantile,probs=c(.05,.95))
          coverage[rep]=sum(sapply(1:length(truth), function(i) truth[i]>=cov.temp[1,i]
                          & truth[i]<=cov.temp[2,i]))/length(truth)
          width[rep]=mean(sapply(1:ncol(cov.temp),function(i) abs(cov.temp[1,i]-cov.temp[2,i])))
        }
        ##printing
        print((sig.holder))
        print(bart.mod$sigest)
        print(width)
        #print(true.tf)
        #print((true.rmse))
        #print(coverage)
        
        big.list[[count]]=list()
        big.list[[count]][[1]]=bart.mod
        big.list[[count]][[2]]=c(tree.vec[i],tf.vec[k],factor.vec[j],mean(sig.holder))
        big.list[[count]][[3]]=mean(prior.sig)
        big.list[[count]][[4]]=mean(true.tf)
        big.list[[count]][[5]]=mean(sel.prop)
        big.list[[count]][[6]]=mean(train.rmse)
        big.list[[count]][[7]]=mean(test.rmse)
        big.list[[count]][[8]]=mean(true.rmse)
        big.list[[count]][[9]]=mean(coverage)
        big.list[[count]][[10]]=mean(width)
        print(count)
        count=count+1
      }
    }
  }
  return(big.list)
}


##3.
bart_sig_post_test=function(tf.exp,gene_list,n.obs,reps,factor.vec,tree.vec,
                       tf.vec,sigest.vec,burn_size,post_size,sig_df=3,
                       sig_quant=.90,thin.size=2,test.list,truth){
  ##Second item in the list has the following order:
  #(tree size, number of tfs, noise factor, posterior sigma estimate)
  ##Need to make more robust to handle different orderings-right now just copies above
  big.list=list()
  count=1
  
  test_factor_idx=sapply(1:length(test.list), function(i) test.list[[i]][[2]])
  test.mat=test.list[[which(is.na(test_factor_idx))]][[1]] ##get test mat
  
  
  for(i in 1:length(tree.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(tf.vec)){
        sig.holder=numeric(reps) ##posterior var vector
        true.tf=numeric(reps) ##vector for correct proportion
        sel.prop=numeric(reps) ##vector for proportion of time selected
        
        train.rmse=numeric(reps) ##in sample RMSE
        test.rmse=numeric(reps) ##vector for RMSE on test set
        true.rmse=numeric(reps) ##distance from truth
        coverage=numeric(reps) ##frequentist coverage
        width=numeric(reps)
        for(rep in 1:reps){
          
          gene.exp=gene_list[[j]][[rep]]
          #print(gene.exp[1:4])
          
          train.exp=tf.exp[1:n.obs,1:tf.vec[k]]
          
          ##Handle n less than p          

          bart.mod = bart(x.train=train.exp,
                          y.train=gene.exp,
                          #x.test=test.mat,
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
      ##Calculations for RMSE
          ##test info
          #train.rmse[reps]=sqrt(sum((bart.mod$yhat.train.mean-gene.exp)^2)/n.obs)
          #loc=which(test_factor_idx==factor.vec[j])
          #y.test=test.list[[loc]][[1]]
          #test.rmse[reps]=sqrt(sum((bart.mod$yhat.test.mean-y.test)^2)/nrow(test.mat))
          true.rmse[rep]=sqrt(sum((bart.mod$yhat.train.mean-truth)^2)/nrow(train.exp))
          cov.temp=apply(bart.mod$yhat.train,2,quantile,probs=c(.05,.95))
          coverage[rep]=sum(sapply(1:length(truth), function(i) truth[i]>=cov.temp[1,i]
                                    & truth[i]<=cov.temp[2,i]))/length(truth)
          width[rep]=mean(sapply(1:ncol(cov.temp),function(i) abs(cov.temp[1,i]-cov.temp[2,i])))
          
        }
        
        
      
        ##printing to check
        print(mean(sig.holder))
        print(bart.mod$sigest)
        print(true.tf)
        
        
        print((true.rmse))
        print((coverage))
        
        big.list[[count]]=list()
        big.list[[count]][[1]]=bart.mod
        big.list[[count]][[2]]=c(tree.vec[i],tf.vec[k],factor.vec[j],mean(sig.holder))
        big.list[[count]][[3]]=NA
        big.list[[count]][[4]]=mean(true.tf)
        big.list[[count]][[5]]=mean(sel.prop)
        big.list[[count]][[6]]=mean(train.rmse)
        big.list[[count]][[7]]=mean(test.rmse)
        big.list[[count]][[8]]=mean(true.rmse)
        big.list[[count]][[9]]=mean(coverage)
        big.list[[count]][[10]]=mean(width)
        
        print(count)
        count=count+1
      }
    }
  }
  return(big.list)
}

##4. 
plot_sig_full=function(bart.list,n_burn,subtitle,y.lims){
  sig=c(bart.list$first.sigma,bart.list$sigma)
  burn.idx=1:n_burn
  burn=sig[burn.idx]
  post=sig[-burn.idx]
  plot(burn,main=paste("Convergence Plot for Sigma\n",subtitle),
       xlab="Sample Number",
       ylab="Sigma",
       xlim=c(0,length(sig)),
       ylim=y.lims,
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