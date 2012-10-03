require(BayesTree)
 setwd("~/Documents/Research/Genomics")
source("bart_fns.R")
source("bart_fns2.R")

##Set-up
set.seed(20)

##TF settings
tf.size1=c(10,15,20,25)
tf.size2=c(10,20,30,40) ##script is set to work with this right now


##Observation settings
n50=50
n100=100
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
n.tree.vec=c(5,10,15,20)
factor.vec1=c(0,.25,.5,1,1.5,2,3.5,5,8)
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12) ##work with this right now 
burn=2500
post=5000

##Prior weight vectors
rep2=c(rep(1,times=2),2:max(tf.size2))
rep5=c(rep(1,times=5),2:max(tf.size2))
rep10=c(rep(1,times=10),2:max(tf.size2))
rep50=c(rep(1,times=50),2:max(tf.size2))


##Get gene data however you need it:
#load("gene_data.RData")
#gene.data=get_gene_data_list(max.n.obs=n300,max.tf.exp=tf.exp.300,
 #                            factor.vec=factor.vec2,tf.beta=tf.beta.1)

##50 observations
output.50.rep2=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.50,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,prior.vec=rep2,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.50.rep5=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.50,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,prior.vec=rep5,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.50.rep10=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.50,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,prior.vec=rep10,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.50.rep50=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.50,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,prior.vec=rep50,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)


#100 Observations
output.100.rep2=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.100,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n100,
                      n.tree.vec=n.tree.vec,prior.vec=rep2,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.100.rep5=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.100,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n100,
                      n.tree.vec=n.tree.vec,prior.vec=rep5,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.100.rep10=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.100,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n100,
                      n.tree.vec=n.tree.vec,prior.vec=rep10,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.100.rep50=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.100,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n100,
                      n.tree.vec=n.tree.vec,prior.vec=rep50,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

##300 Observations
output.300.rep2=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.300,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n300,
                      n.tree.vec=n.tree.vec,prior.vec=rep2,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.300.rep5=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.300,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n300,
                      n.tree.vec=n.tree.vec,prior.vec=rep5,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.300.rep10=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.300,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n300,
                      n.tree.vec=n.tree.vec,prior.vec=rep10,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.300.rep50=bart_loop_prior(gene.list=gene.data,tf.exp=tf.exp.300,tf.beta=tf.beta.1,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n300,
                      n.tree.vec=n.tree.vec,prior.vec=rep50,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)


#Save Files
save(output.50.rep2,file="rep2_50.R")
save(output.50.rep5,file="rep5_50.R")
save(output.50.rep10,file="rep10_50.R")
save(output.50.rep50,file="rep50_50.R")

save(output.100.rep2,file="rep2_100.R")
save(output.100.rep5,file="rep5_100.R")
save(output.100.rep10,file="rep10_100.R")
save(output.100.rep50,file="rep50_100.R")

save(output.300.rep2,file="rep2_300.R")
save(output.300.rep5,file="rep5_300.R")
save(output.300.rep10,file="rep10_300.R")
save(output.300.rep50,file="rep50_300.R")

##Plots
bart_plot_both_prior(output.50.rep2,n.tree.vec,factor.vec2,tf.size2,50,0,2)
bart_plot_both_prior(output.50.rep5,n.tree.vec,factor.vec2,tf.size2,50,0,5)
bart_plot_both_prior(output.50.rep10,n.tree.vec,factor.vec2,tf.size2,50,0,10)
bart_plot_both_prior(output.50.rep50,n.tree.vec,factor.vec2,tf.size2,50,0,50)

bart_plot_both_prior(output.100.rep2,n.tree.vec,factor.vec2,tf.size2,100,0,2)
bart_plot_both_prior(output.100.rep5,n.tree.vec,factor.vec2,tf.size2,100,0,5)
bart_plot_both_prior(output.100.rep10,n.tree.vec,factor.vec2,tf.size2,100,0,10)
bart_plot_both_prior(output.100.rep50,n.tree.vec,factor.vec2,tf.size2,100,0,50)

bart_plot_both_prior(output.300.rep2,n.tree.vec,factor.vec2,tf.size2,300,0,2)
bart_plot_both_prior(output.300.rep5,n.tree.vec,factor.vec2,tf.size2,300,0,5)
bart_plot_both_prior(output.300.rep10,n.tree.vec,factor.vec2,tf.size2,300,0,10)
bart_plot_both_prior(output.300.rep50,n.tree.vec,factor.vec2,tf.size2,300,0,50)

##Summaries
##50 obs
summary(output.50.1$"True_TF_Prop");sd(output.50.1$"True_TF_Prop")
sum(output.50.1$"True_Most_Common?")/nrow(output.50.1)

summary(output.50.rep2$"True_TF_Prop");sd(output.50.rep2$"True_TF_Prop")
sum(output.50.rep2$"True_Most_Common?")/nrow(output.50.rep2)

summary(output.50.rep5$"True_TF_Prop");sd(output.50.rep5$"True_TF_Prop")
sum(output.50.rep5$"True_Most_Common?")/nrow(output.50.rep5)

summary(output.50.rep10$"True_TF_Prop");sd(output.50.rep10$"True_TF_Prop")
sum(output.50.rep10$"True_Most_Common?")/nrow(output.50.rep10)

summary(output.50.rep50$"True_TF_Prop");sd(output.50.rep50$"True_TF_Prop")
sum(output.50.rep50$"True_Most_Common?")/nrow(output.50.rep50)

##100 Obs
summary(output.100.1$"True_TF_Prop");sd(output.100.1$"True_TF_Prop")
sum(output.100.1$"True_Most_Common?")/nrow(output.100.1)

summary(output.100.rep2$"True_TF_Prop");sd(output.100.rep2$"True_TF_Prop")
sum(output.100.rep2$"True_Most_Common?")/nrow(output.100.rep2)

summary(output.100.rep5$"True_TF_Prop");sd(output.100.rep5$"True_TF_Prop")
sum(output.100.rep5$"True_Most_Common?")/nrow(output.100.rep5)

summary(output.100.rep10$"True_TF_Prop");sd(output.100.rep10$"True_TF_Prop")
sum(output.100.rep10$"True_Most_Common?")/nrow(output.100.rep10)
                                              
summary(output.100.rep50$"True_TF_Prop");sd(output.100.rep50$"True_TF_Prop")
sum(output.100.rep50$"True_Most_Common?")/nrow(output.100.rep50)                                         

##300 Obs
summary(output.300.1$"True_TF_Prop");sd(output.300.1$"True_TF_Prop")
sum(output.300.1$"True_Most_Common?")/nrow(output.300.1)

summary(output.300.rep2$"True_TF_Prop");sd(output.300.rep2$"True_TF_Prop")
sum(output.300.rep2$"True_Most_Common?")/nrow(output.300.rep2)

summary(output.300.rep5$"True_TF_Prop");sd(output.300.rep5$"True_TF_Prop")
sum(output.300.rep5$"True_Most_Common?")/nrow(output.300.rep5)

summary(output.300.rep10$"True_TF_Prop");sd(output.300.rep10$"True_TF_Prop")
sum(output.300.rep10$"True_Most_Common?")/nrow(output.300.rep10)

summary(output.300.rep50$"True_TF_Prop");sd(output.300.rep50$"True_TF_Prop")
sum(output.300.rep50$"True_Most_Common?")/nrow(output.300.rep50)


#Functions
bart_loop_prior=function(gene.list,tf.exp,tf.beta,factor.vec,tf.size.vec,n.obs,
                         n.tree.vec,prior.vec,burn_size,
                         post_size,thin.size=1,sigest=NA){
  output=data.frame()
  count=1
  gene_factor_idx=sapply(1:length(gene.list), function(i) gene.list[[i]][[2]]) ##returns indices of noise level data set was generated by
  for(i in 1: length(tf.size.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(n.tree.vec)){
        loc=which(gene_factor_idx==factor.vec[j])
        gene.exp=gene.list[[loc]][[1]][1:n.obs] ##returns first n.obs elements
        print(c(loc,gene.exp[1:4]))
        
        ##Set up number of repititions and repeats in data frame
        rep.nums=sapply(1:max(tf.size.vec),function(x) length(which(prior.vec==x)))
        #print(rep.nums)
        reps=rep.nums[1:tf.size.vec[i]]
        #print(reps)
        train.exp.temp=tf.exp[1:n.obs,1:tf.size.vec[i]]
        train.exp=train.exp.temp[,rep(1:ncol(train.exp.temp),reps)]
        #
        
        ##For case of N-use sample SD-else if will use whatever is pre-specified
        if(n.obs<=ncol(train.exp) & is.na(sigest)) bart.sig=sd(gene.exp)
        else bart.sig=NA ##I think this should be sigest-check later
        print(dim(train.exp))
        print(bart.sig)
        #print(train.exp[1:2,])
        
        ##Actual BART Runs
        
        bart.mod = bart(x.train=train.exp,
              y.train=gene.exp,
                ntree=n.tree.vec[k],
                sigest=bart.sig,
                nskip=burn_size,
                ndpost=post_size,
                keepevery=thin.size,
                verbose=F) ##Suppress printing
        
        ##Start creating output               
        ##Basic Storage
        output[count,1]=tf.size.vec[i] ##First location is number of transcription factors
        output[count,2]=factor.vec[j] ##Amount of noise
        output[count,3]=n.tree.vec[k]
        
        last.spot=max(which(prior.vec==tf.size.vec[i])) ##computes last spot in 
                      ##prior vector that will be needed for sum and prop calcs
        
        print(last.spot)
        #print(sum_calc(bart.mod))
        sums=sum_calc_prior(bart.mod,prior.vec[1:last.spot])
        #print(sums)
        #print(sum(sums))
        #props=prop_calc(bart.mod)
        props=prop_calc_prior(bart.mod,prior.vec[1:last.spot])              
        #print(props)
        sorted.sums=sort(sums,decreasing=T)
        sorted.props=sort(props,decreasing=T)    
        #print(sorted.sums)
        
        output[count,4]=which.max(sums)==1  ##Does our appear most?               
        output[count,5]=sorted.props[1] ##Top three TFs
        output[count,6]=which(props==sorted.props[1])[1] ##How to deal with ties?
        output[count,7]=sorted.props[2]
        output[count,8]=which(props==sorted.props[2])[1]
        output[count,9]=sorted.props[3]
        output[count,10]=which(props==sorted.props[3])[1]
        output[count,11]=sums[1]
        output[count,12]=sum(sums) ##Total splits 
        output[count,13]=round(bart.mod$sigest,3)
        output[count,14]=props[1]
        print(count)                 
        count=count+1                    
      }               
    }
  }  
  colnames(output)=c("Num_TFs","Noise/Signal","Num_Trees","True_Most_Common?","1st","Name-1st",
                  "2nd","Name-2nd","3rd","Name-3rd","True_TF_Sum","Tot_Splits","Sigma_Estimate","True_TF_Prop")
  return(output)
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

#4.
##Plotting function for holding number of trees constant.
##Takes output from bart_loop()
bart_plot_const_tree_prior=function(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent,prior){
  ##Initial Plot Set-up
    par(mgp=c(1.8,.5,0))
  col.vec=c("orange","green","blue","red")
  plot(1,type="n",xlim=c(min(ntree.vec),max(ntree.vec)+4),ylim=c(min(factor.vec),max(factor.vec)+1),
       xlab="Number of Trees",ylab="Noise/Signal Ratio",
       main=paste("Break Levels for BART Models\n",num.obs,"Observations,", paste(prior, "x",sep=""), "Weight"))
  
  for(i in 1:length(ntree.vec)){
    for(j in 1:length(tfsize.vec)){
      temp=subset(output,subset=Num_Trees==ntree.vec[i] & Num_TFs==tfsize.vec[j])
      noise_count=1
      repeat{
      if(temp$"True_Most_Common?"[noise_count]!=T || 
        (temp$"1st"[noise_count]-temp$"2nd"[noise_count])/temp$"2nd"[noise_count]<percent) {
          noise_count=noise_count
          break
      } 
            if(noise_count==nrow(temp)) break
            noise_count=noise_count+1
            
      }
  noise_count    
   points(jitter(temp$"Num_Trees"[noise_count],factor=.3),jitter(temp$"Noise/Signal"[noise_count],factor=2.5),pch=(14+j),col=col.vec[j]) 
    }
  }
      leg.names=paste(tfsize.vec,"TFs")
     legend("topright",legend=leg.names,pch=(14+1):(14+length(tfsize.vec)),col=col.vec,cex=.8)
}

#5.
##Plotting function for holding number of TFs constant.
##Takes output from bart_loop()
bart_plot_const_tf_prior=function(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent, prior){
  ##Initial Plot Setup
    par(mgp=c(1.8,.5,0))
  col.vec=c("orange","green","blue","red")
  plot(1,type="n",ylim=c(min(factor.vec),max(factor.vec)+1),xlim=c(min(tfsize.vec),max(tfsize.vec)+8),
       xlab="Number of TFs",ylab="Noise/Signal Ratio",
       main=paste("Break Levels for BART Models\n",num.obs,"Observations,", paste(prior, "x",sep=""), "Weight"),
                  "xaxt"="n")
  axis(side=1,at=tfsize.vec)
  for(i in 1:length(tfsize.vec)){
    for(j in 1:length(ntree.vec)){
      temp=subset(output,subset=Num_Trees==ntree.vec[j] & Num_TFs==tfsize.vec[i])
      noise_count=1
      repeat{
      if(temp$"True_Most_Common?"[noise_count]!=T || 
        (temp$"1st"[noise_count]-temp$"2nd"[noise_count])/temp$"2nd"[noise_count]<percent) {
          noise_count=noise_count
          break
      } 
            if(noise_count==nrow(temp)) break
            noise_count=noise_count+1
            
      }    
   points(jitter(temp$"Num_TFs"[noise_count],factor=.3),jitter(temp$"Noise/Signal"[noise_count],factor=2.5),pch=(14+j),col=col.vec[j]) 
    }
  }
      leg.names=paste(ntree.vec,"Trees")
     legend("topright",legend=leg.names,pch=(14+1):(14+length(tfsize.vec)),col=col.vec,cex=.8)
}


##6.
##Plots both views-holding TF and tree constant 
bart_plot_both_prior=function(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent,prior){
  par(mfrow=c(1,2))
  bart_plot_const_tree_prior(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent,prior)
    bart_plot_const_tf_prior(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent,prior)
}



####Test Stuff
a=c(1,1,2,2,3,3)
max(which(a==3))
v=cbind(bart.mod$varcount[,1],bart.mod$varcount)
c=apply(v,2,sum)
length(c)
c
tapply(c,rep2,sum)
length(rep2)
aggregate(v,by=as.list(rep2),sum)
ag
colnames(v)

t=cbind(tf.exp.30[,1],tf.exp.30[,1],tf.exp.30)
g=tf.exp.30%*%tf.beta.1

sum(sum_calc(bart.mod))
prop_calc_prior(bart.mod,rep3[1:12])
rep3=c(rep(1,times=3),2:40)

  bart.mod = bart(x.train=t[,1:12],
              y.train=g,
                ntree=5,
                sigest=sd(g),
                nskip=300,
                ndpost=1000,
                keepevery=2,
                verbose=F) ##Suppress printing


