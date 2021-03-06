##Script for generating test data
require(BayesTree)
setwd("~/Documents/Research/Genomics/BART_2/Code_Objects")
source("bart_fns.R")
source("bart_fns2_2.R")

##Set-up
set.seed(20)

##TF settings
tf.size1=c(10,15,20,25)
tf.size2=c(10,20,30,40) ##script is set to work with this right now


##Observation settings
n30=30
n50=50
n100=100
n200=200
n300=300

mean.tf=0 #1.4862e-05 ##generate the X matrix which will be fixed- sl
sd.tf=1 #.4778817-not using just to keep larger numbers at play

##Generate Design Matrices

tf.exp.300=sapply(rep(n300,max(tf.size2)),rnorm,mean=mean.tf,sd=sd.tf) #gives full matrix
tf.exp.30=tf.exp.300[1:n30,]
tf.exp.50=tf.exp.300[1:n50,]
tf.exp.100=tf.exp.300[1:n100,]
tf.exp.200=tf.exp.300[1:n200,]

##Beta settings
tf.beta.1=c(1,rep(0,times=max(tf.size2-1)))

##signal-300
signal.300.1=sum(abs(tf.exp.300%*%tf.beta.1))/n300
signal.300.1

##Function params
n.tree.vec=c(5,10,15,20)
factor.vec1=c(0,.25,.5,1,1.5,2,3.5,5,8)
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12) ##work with this right now 
burn=2500
post=5000



##Function returns a list of lists. The first element in each list is a vector for a
##specific noise level. The second is the corresponding level of noise.
##The last element contains the test matrix and can be called by which(is.na()) through
##an sapply loop.
##Also need to put in the signal from the matrix used to generate the noise.
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

##Inputs
signal.300.1=sum(abs(tf.exp.300%*%tf.beta.1))/n300
n.obs=250
n.cols=40
tf.beta=tf.beta.1=c(1,rep(0,times=max(n.cols-1)))
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12)
#output
test.list=gen_test(250,40,tf.beta,factor.vec2,signal.300.1)

sapply(1:length(test.list), function(i) test.list[[i]][[2]])
dim(test.list[[18]][[1]])
save(test.list,file="test_set.RData")



#3.
##Main function that calculates different statistics for BART simulations-single gene
##re-written to use same y-vectors (only different for noise)
##now returns list with 2 elements. 1st is usual BART output from bart_loop function
##and second is information about predictive performance

##Rewritten to compute sqrt(||XB-yhat||^2/n)
bart_loop_test=function(gene.list,tf.exp,tf.beta,factor.vec,tf.size.vec,n.obs,
                   n.tree.vec,burn_size,post_size,thin.size=1,sigest=NA,base=.9,
                        test.list,signal){
  out.list=list()
  output=data.frame()
  test.output=data.frame()
  count=1
  gene_factor_idx=sapply(1:length(gene.list), function(i) gene.list[[i]][[2]]) ##returns indices of noise level data set was generated by
  ##setup for test
  #test_factor_idx=sapply(1:length(test.list), function(i) test.list[[i]][[2]])
  #test_mat_idx=which(is.na(test_factor_idx))
  #test.mat.full=test.list[[test_mat_idx]][[1]]
  #test.n.obs=nrow(test.mat)
  truth=tf.exp%*%tf.beta
  
  for(i in 1: length(tf.size.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(n.tree.vec)){
        loc=which(gene_factor_idx==factor.vec[j])
        gene.exp=gene.list[[loc]][[1]][1:n.obs] ##returns first n.obs elements
        print(c(loc,gene.exp[1:4]))
        ##of data vector for corresponding noise
        train.exp=tf.exp[1:n.obs,1:tf.size.vec[i]] ##adjust size
        #test.mat=test.mat.full[,1:tf.size.vec[i]] ##adjust size to match num cols
        ##For case of N-use sample SD-else if will use whatever is pre-specified
        if(n.obs<=ncol(train.exp) & is.na(sigest)) bart.sig=sd(gene.exp)
        else bart.sig=NA
        print(dim(train.exp))
        print(bart.sig)
        bart.mod = bart(x.train=train.exp,
                        y.train=gene.exp,
                        #x.test=test.mat,
                        ntree=n.tree.vec[k],
                        sigest=bart.sig,
                        nskip=burn_size,
                        ndpost=post_size,
                        keepevery=thin.size,
                        base=base,
                        verbose=F) ##Suppress printing
        
        output[count,1]=tf.size.vec[i] ##First location is number of transcription factors
        output[count,2]=factor.vec[j] ##Amount of noise
        output[count,3]=n.tree.vec[k]
        
        sums=sum_calc(bart.mod)
        props=prop_calc(bart.mod)
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
        output[count,13]=round(bart.mod$sigest,3)
        output[count,14]=props[1]
        
        test.output[count,1]=tf.size.vec[i] ##First location is number of transcription factors
        test.output[count,2]=factor.vec[j] ##Amount of noise
        test.output[count,3]=n.tree.vec[k]
        test.output[count,4]=signal*factor.vec[j] #sigma
        #Compute true RMSE
        test.output[count,6]=sqrt(sum((bart.mod$yhat.train.mean-truth)^2)/n.obs)
        #Compute test-set RMSE
        #test.gene.exp=test.list[[which(test_factor_idx==factor.vec[j])]][[1]]
        #test.output[count,5]=sqrt(sum((bart.mod$yhat.test.mean-test.gene.exp)^2)/test.n.obs)
        
#         test.output[[count]]=list()
#         test.output[[count]][["train"]]=bart.mod$yhat.train.mean
#         test.output[[count]][["test"]]=bart.mod$yhat.test.mean
#         test.output[[count]][["factor"]]=factor.vec[j]
#         test.output[[count]]
#         test.output[[count]][["y_vec"]]=
        print(count)                 
        count=count+1                    
      }               
    }
  }  
  colnames(output)=c("Num_TFs","Noise/Signal","Num_Trees","True_Most_Common?","1st","Name-1st",
                     "2nd","Name-2nd","3rd","Name-3rd","True_TF_Sum","Tot_Splits","Sigma_Estimate","True_TF_Prop")
  
  colnames(test.output)=c("Num_TFs","Noise/Signal","Num_Trees","True_Sigma","True_RMSE")
  out.list[[1]]=output
  out.list[[2]]=test.output
  return(out.list)
}


test.50.1=bart_loop_test(gene.list=gene.data,tf.exp=tf.exp.50,tf.beta=tf.beta.1,
                      factor.vec=factor.vec2,tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA,test.list=test.list,
                         signal=signal.300.1)


##Junk

test.mat=test.list[[length(test.list)]][[1]]
gene.vec=gene.data[[1]][[1]]
gene.vec
bm=bart(x.train=tf.exp.300,y.train=gene.vec,x.test=test.mat,ntree=100,ndpost=500,nskip=500)
length(bm$yhat.train.mean)
sum((gene.vec-bm$yhat.train.mean)^2)/300
sum((bm$yhat.test.mean-d)^2)/250
x=sum(bm$yhat.test.mean-d)
y=x^2
sum(y)/250
d=test.list[[1]][[1]]
p=(lm(d~test.mat))
summary(p)
sqrt(sum((p$fitted.values-d)^2)/209)
p$fitted.values
anova(p)
sqrt(sum((d-mean(d))^2)/250)
sum(mean(test.list[[1]][[1]])-test.list[[1]][[1]])^2/250
length(bm$yhat.train.mean)
test.list[[1]][[2]]

sum(abs(test.mat[,1]))/250
