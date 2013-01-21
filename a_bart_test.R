setwd("~/Research_Genomics/bart_gene") ##work
setwd("~/Documents/Research/Genomics/bart_gene/") ##home
source("real_functions.R")
library(BayesTree)

##work
priors=read.table("~/Research_Genomics/CHIP.priorprobs.39.txt",header=T) ##work
gene.exp=read.table("~/Research_Genomics/expression.genes.txt",header=T) ##work
tf.exp=read.table("~/Research_Genomics/expression.tfs.39.txt",header=T) ##work
##Home
priors=read.table("~/Documents/Research/Genomics/Real_Data/CHIP.priorprobs.39.txt",header=T) ##home
gene.exp=read.table("~/Documents/Research/Genomics/Real_Data/expression.genes.txt",header=T) ##home
tf.exp=read.table("~/Documents/Research/Genomics/Real_Data/expression.tfs.39.txt",header=T) ##home

##put in setup
priorMat=as.matrix(priors[,3:ncol(priors)])
rownames(priorMat)=priors[,1] ; colnames(priorMat)=colnames(priors)[3:ncol(priors)]
out=setup() ##5 objects
priorWeights=out[[1]]#PriorWeights: rows are genes, cols are TFs
gene.train=out[[4]];gene.test=out[[5]]##gene.train and gene.test rows are obs and cols are genes
tf.train=out[[2]];tf.test=out[[3]]##TF train and TF test: rows are obs and cols are TFs
geneNames=as.character(gene.exp[,2]) ##gene names

##BART setup stuff
setwd("C:/Users/jbleich/workspace/CGMBART_GPL")
directory_where_code_is = "C:\\Users\\jbleich\\workspace\\CGMBART_GPL" ##hack because of bart package
source("r_scripts/bart_package.R")

##for grant
nrows=250
ncols=20
X=matrix(rnorm(nrows*ncols,0,1),nrows,ncols)
y=3*X[,1]+2*X[,2]+rnorm(nrows,0,1)
X_test=matrix(rnorm(nrows*ncols,0,1),nrows,ncols)
y_test=3*X_test[,1]+2*X_test[,2]+rnorm(nrows,0,1)
data=data.frame(cbind(X,y))
test=data.frame(X_test,y_test)

prior_vec=c(100,100,rep(1,times=ncols-2))

bart_machine=build_bart_machine(training_data=data,num_trees=200,num_burn_in=2000,num_iterations_after_burn_in=5000,run_in_sample=T)
sqrt(bart_machine$L2_err_train/nrows)
preds=bart_predict(bart_machine=bart_machine,new_data=as.data.frame(X_test))
sqrt(sum((preds-y_test)^2)/nrows)

bart_machine_prior=build_bart_machine(training_data=data,num_trees=200,num_burn_in=2000,num_iterations_after_burn_in=5000,run_in_sample=T,cov_prior_vec=prior_vec)
sqrt(bart_machine_prior$L2_err_train/nrows)
preds_prior=bart_predict(bart_machine=bart_machine_prior,new_data=as.data.frame(X_test))
sqrt(sum((preds_prior-y_test)^2)/nrows)

rf=randomForest(x=X,y=y)
train_rf=predict(object=rf,newdata=X)
sqrt(sum((y-train_rf)^2)/nrows)
preds_rf=predict(object=rf,newdata=X_test)
sqrt(sum((y_test-preds_rf)^2)/nrows)
##

#########old Dec. 2012#################
pvec=priors[2,]
pvec[which(pvec==0)]=.001
pvec=pvec[3:length(pvec)]
pvecs=pvec[sample(1:39,39,replace=F)]
y=getGeneResponse(geneNames[2])
bart.training=data.frame(tf.train,y)
bartmod=build_bart_machine(training_data=bart.training,num_trees=100,num_burn_in=2000,num_iterations_after_burn_in=2000,run_in_sample=T)
bartmod$L2_err_train
bartmod=build_bart_machine(training_data=bart.training,num_trees=100,num_burn_in=2000,num_iterations_after_burn_in=2000)
bartmod2=build_bart_machine(training_data=bart.training,num_trees=200,num_burn_in=2000,num_iterations_after_burn_in=2000)
counts=get_var_counts_over_chain(bartmod)
props=apply(counts,2,sum)/sum(apply(counts,2,sum))
names(props)=names(pvec)
props=round(props,3)
sort(pvec,decreasing=T)
sort(props,decreasing=T)
plot(props,pvec)
counts2=get_var_counts_over_chain(bartmod2)
props2=apply(counts2,2,sum)/sum(apply(counts2,2,sum))
names(props2)=names(pvec)
props2=round(props2,3)
sort(props,decreasing=T)
bayes=bart(x.train=tf.train,y.train=y,ntree=200,nskip=2000,ndpost=2000)
bayes2=bart(x.train=tf.train,y.train=y,ntree=200,nskip=2000,ndpost=2000)
pbayes=prop_calc(bayes)
names(pbayes)=names(pvec)
pbayes2=prop_calc(bayes2)
names(pbayes2)=names(pvec)
sort(pbayes,decreasing=T)
sort(props,decreasing=T)
plot(props,pbayes)
plot(pbayes2,pbayes)
plot(props,props2)
props
props2
pbayes
pbayes2
bp=bart_predict_for_test_data(bart_machine=bartmod,test_data=bart.training)
mean(bp$y_hat)
bp2=bart_predict_for_test_data(bart_machine=bartmod2,test_data=bart.training)
mean(bp2$y_hat)
mean(bayes$yhat.train.mean)
mean(bayes2$yhat.train.mean)
plot(bp$y_hat,bp2$y_hat)
plot(bayes$yhat.train.mean,bayes2$yhat.train.mean)

 nrows=250
ncols=20
x.train=matrix(rnorm(nrows*ncols,0,1),nrow=nrows,ncol=ncols)
x.test=matrix(rnorm(nrows*ncols,0,1),nrow=nrows,ncol=ncols)
y=rnorm(250)
data=data.frame(x.train,y)
bart_machine=build_bart_machine(training_data=data,num_trees=20,num_burn_in=1000,num_iterations_after_burn_in=2000)
num_cores=1
var_num=1

sig=get_variable_significance(bart_machine,data=data,var_num=1,num_iter=100,print_histogram=T)
sig$sse_vec
sig$real_sse



beta=c(1,2,rep(0,18))
y=x.train%*%beta+rnorm(250,0,1)
y.test=x.test%*%beta+rnorm(250,0,1)
prior=c(1,1,rep(.001,18))
data2=data.frame(x.train,y)
bayes=bart(x.train=x.train,y.train=as.numeric(y),x.test=x.test,ntree=200,nskip=1000,ndpost=2000)
sum((y-bayes$yhat.train.mean)^2)
sum((y.test-bayes$yhat.test.mean)^2)

bart_machine=build_bart_machine(training_data=data2,num_trees=50,num_burn_in=1000,num_iterations_after_burn_in=2000,run_in_sample=T,cov_prior_vec=prior)
bart_machine$L2_err_train
apply(get_var_counts_over_chain(bart_machine),2,sum)
plot(y,bart_machine$y_hat)
test.data=data.frame(x.test,y.test);colnames(test.data)[21]="y"
pred=bart_predict_for_test_data(bart_machine=bart_machine,test_data=test.data)
pred$L2_err
rf=randomForest(x=x.train,y=y,ntree=500)
prf=predict(rf,x.test)
sum((y.test-prf)^2)

sig1=get_variable_significance(bart_machine,data=data2,var_num=1,num_iter=100,print_histogram=T)
sig3=get_variable_significance(bart_machine,data=data2,var_num=3,num_iter=100,print_histogram=T)
sig1$real_sse-mean(sig1$sse_vec)

fullmod = lm(y~., data = data2)
redmod = lm(y~X1+X2, data = data2)
anova(fullmod, redmod)

sig3$real_sse-mean(sig3$sse_vec)
sig3$real