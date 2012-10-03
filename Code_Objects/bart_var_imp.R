  ##Variable Importance Plots
  
  tb=tf.beta.1
  mat=tf.exp.300[1:10,]
  test.mat=sapply(rep(100,40),rnorm,mean=0,sd=1) 
  
  shuffle=sapply(1:2,function(i) sample.int(nrow(test.mat)) ,simplify=T)
  shuff.mat1=test.mat
  shuff.mat2=test.mat
  shuff.mat1[,1]=shuff.mat1[shuffle[,1],1]  
  shuff.mat2[,2]=shuff.mat2[shuffle[,2],2]                
  y.vec=as.numeric(mat%*%tb+rnorm(nrow(mat),0,1))
  y.test=as.numeric(test.mat%*%tb+rnorm(nrow(test.mat),0,1))
  
  set.seed(15)
  bart.mod0=bart(x.train=mat,y.train=y.vec,x.test=test.mat,
                 ntree=20,ndpost=1000,nskip=1000,sigest=sd(y.vec))                
                  
  set.seed(15)
  bart.mod1=bart(x.train=mat,y.train=y.vec,x.test=shuff.mat1,
                 ntree=20,ndpost=1000,nskip=1000,sigest=sd(y.vec))                
                  
  set.seed(15)
  bart.mod2=bart(x.train=mat,y.train=y.vec,x.test=shuff.mat2,
                  ntree=20,ndpost=1000,nskip=1000,sigest=sd(y.vec))                
                  
  yh0=bart.mod0$yhat.train.mean                            
  yh1=bart.mod1$yhat.train.mean
  yh2=bart.mod2$yhat.train.mean
  
  t0=bart.mod0$yhat.test.mean          
  t1=bart.mod1$yhat.test.mean
  t2=bart.mod2$yhat.test.mean
  
  r0=sqrt(sum((t0-y.test)^2)/nrow(test.mat))
  r1=sqrt(sum((t1-y.test)^2)/nrow(test.mat))
  r2=sqrt(sum((t2-y.test)^2)/nrow(test.mat))
  
  
  dis0=sqrt(sum((yh0-mat%*%tb)^2)/nrow(mat))
  dis1=sqrt(sum((yh1-mat%*%tb)^2)/nrow(mat))
  dis2=sqrt(sum((yh2-mat%*%tb)^2)/nrow(mat))
  
  dis0
  dis1
  dis2
  
  r0
  r1
  r2                
                  
                  
      
                  
                
                
                
                
                