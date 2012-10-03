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
tf.beta.2=c(2,rep(0,times=max(tf.size2-1)))

#signal setting-30
signal.30.1=sum(abs(tf.exp.30%*%tf.beta.1))/n30
signal.30.1
signal.30.2=sum(abs(tf.exp.30%*%tf.beta.2))/n30
signal.30.2
#signal setting-50
signal.50.1=sum(abs(tf.exp.50%*%tf.beta.1))/n50
signal.50.1
signal.50.2=sum(abs(tf.exp.50%*%tf.beta.2))/n50
signal.50.2
##signal-100
signal.100.1=sum(abs(tf.exp.100%*%tf.beta.1))/n100
signal.100.1
##signal-300
signal.200.1=sum(abs(tf.exp.200%*%tf.beta.1))/n200
signal.200.1
##signal-300
signal.300.1=sum(abs(tf.exp.300%*%tf.beta.1))/n300
signal.300.1

##Function params
n.tree.vec=c(5,10,15,20)
factor.vec1=c(0,.25,.5,1,1.5,2,3.5,5,8)
factor.vec2=c(.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12) ##work with this right now 
burn=2500
post=5000



##Get gene data
#save(gene.data,file="gene_data.RData")
gene.data=get_gene_data_list(max.n.obs=n300,max.tf.exp=tf.exp.300,
                             factor.vec=factor.vec2,tf.beta=tf.beta.1)

##Generate Output Data Frames
output.30.1=bart_loop(gene.list=gene.data,tf.exp=tf.exp.30,tf.beta=tf.beta.1,
                      factor.vec=factor.vec2,tf.size.vec=tf.size2,n.obs=n30,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.50.1=bart_loop(gene.list=gene.data,tf.exp=tf.exp.50,tf.beta=tf.beta.1,
                      factor.vec=factor.vec2,tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.100.1=bart_loop(gene.list=gene.data,tf.exp=tf.exp.100,tf.beta=tf.beta.1,
                       factor.vec=factor.vec2,tf.size.vec=tf.size2,n.obs=n100,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                       post_size=post,thin.size=2,sigest=NA)

output.200.1=bart_loop(gene.list=gene.data,tf.exp=tf.exp.200,tf.beta=tf.beta.1,
                       factor.vec=factor.vec2,tf.size.vec=tf.size2,n.obs=n200,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                       post_size=post,thin.size=2,sigest=NA)

output.300.1=bart_loop(gene.list=gene.data,tf.exp=tf.exp.300,tf.beta=tf.beta.1,
                       factor.vec=factor.vec2,tf.size.vec=tf.size2,n.obs=n300,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                       post_size=post,thin.size=2,sigest=NA)


##Second set for beta=2
##Generate Output Data Frames
output.30.2=bart_loop(tf.exp=tf.exp.30,tf.beta=tf.beta.2,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n30,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.50.2=bart_loop(tf.exp=tf.exp.50,tf.beta=tf.beta.2,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n50,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                      post_size=post,thin.size=2,sigest=NA)

output.100.2=bart_loop(tf.exp=tf.exp.100,tf.beta=tf.beta.2,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n100,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                       post_size=post,thin.size=2,sigest=NA)

output.200.2=bart_loop(tf.exp=tf.exp.200,tf.beta=tf.beta.2,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n200,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                       post_size=post,thin.size=2,sigest=NA)

output.300.2=bart_loop(tf.exp=tf.exp.300,tf.beta=tf.beta.2,factor.vec=factor.vec2,
                      tf.size.vec=tf.size2,n.obs=n300,
                      n.tree.vec=n.tree.vec,burn_size=burn,
                       post_size=post,thin.size=2,sigest=NA)



##Second Set to Save-files with 10-40 TFs and larger set of factors
save(output.30.1,file="bart301_v1_2.R")
save(output.50.1,file="bart501_v1_2.R")
save(output.100.1,file="bart1001_v1_2.R")
save(output.200.1,file="bart2001_v1_2.R")
save(output.300.1,file="bart3001_v1_2.R")


bart_plot_both(output.30.1,n.tree.vec,factor.vec2,tf.size2,30,0)
bart_plot_both(output.50.1,n.tree.vec,factor.vec2,tf.size2,50,0)
bart_plot_both(output.100.1,n.tree.vec,factor.vec2,tf.size2,100,0)
bart_plot_both(output.200.1,n.tree.vec,factor.vec2,tf.size2,200,0)
bart_plot_both(output.300.1,n.tree.vec,factor.vec2,tf.size2,300,0)
output.300.1

##Old Code########################################################
#Output for 7_26 report
#save(output.30.1,file="output301.R")
#save(output.30.2,file="output302.R")
#save(output.50.1,file="output501.R")
#save(output.50.2,file="output502.R")

##Old Plotting Commands-for 7_26 report
#bart_loop_plot(5,output.30.1,tf.size)
#bart_loop_plot(10,output.30.1,tf.size)
#bart_loop_plot(15,output.30.1,tf.size)
#bart_loop_plot(20,output.30.1,tf.size)

# bart_loop_plot(5,output.30.2,tf.size)
# bart_loop_plot(10,output.30.2,tf.size)
# bart_loop_plot(15,output.30.2,tf.size)
# bart_loop_plot(20,output.30.2,tf.size)
# 
# bart_loop_plot(5,output.50.1,tf.size)
# bart_loop_plot(10,output.50.1,tf.size)
# bart_loop_plot(15,output.50.1,tf.size)
# bart_loop_plot(20,output.50.1,tf.size)
# 
# bart_loop_plot(5,output.50.2,tf.size)
# bart_loop_plot(10,output.50.2,tf.size)
# bart_loop_plot(15,output.50.2,tf.size)
# bart_loop_plot(20,output.50.2,tf.size)
# 
# bart_loop_plot(5,output.300.1,tf.size)
# bart_loop_plot(10,output.300.1,tf.size)
# bart_loop_plot(15,output.300.1,tf.size)
# bart_loop_plot(20,output.300.1,tf.size)
    