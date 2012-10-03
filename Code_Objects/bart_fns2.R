##Functions
##BART Functions 


#1.
##calculate variable sum
sum_calc=function(bart.obj){apply(bart.obj$varcount,2,sum)}

#2.
##calculate variable proportions 
prop_calc=function(bart.obj){
  tot=sum(apply(bart.obj$varcount,2,sum))
  props=apply(bart.obj$varcount,2,sum)/tot
  round(props,3)
  }

#3.
##Main function that calculates different statistics for BART simulations-single gene
bart_loop=function(tf.exp,tf.beta,factor.vec,tf.size.vec,n.obs,n.tree.vec,burn_size,post_size,thin.size=1,sigest=NA){
  output=data.frame()
  count=1
  for(i in 1: length(tf.size.vec)){
    for(j in 1:length(factor.vec)){
      for(k in 1:length(n.tree.vec)){
        signal.temp=sum(abs(tf.exp%*%tf.beta))/n.obs
        sigma=signal.temp*factor.vec[j]
        gene.exp=as.numeric(tf.exp%*%tf.beta+rnorm(n.obs,mean=0,sd=sigma))
        train.exp=tf.exp[1:n.obs,1:tf.size.vec[i]]
        ##For case of N-use sample SD-else if will use whatever is pre-specified
        if(n.obs<=ncol(train.exp) & is.na(sigest)) bart.sig=sd(gene.exp)
        else bart.sig=NA
        print(dim(train.exp))
        print(bart.sig)
        bart.mod = bart(x.train=train.exp,
              y.train=gene.exp,
                ntree=n.tree.vec[k],
                sigest=bart.sig,
                nskip=burn_size,
                ndpost=post_size,
                keepevery=thin.size,
                base=.65,
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
        print(count)                 
        count=count+1                    
      }               
    }
  }  
  colnames(output)=c("Num_TFs","Noise/Signal","Num_Trees","True_Most_Common?","1st","Name-1st",
                  "2nd","Name-2nd","3rd","Name-3rd","True_TF_Sum","Tot_Splits","Sigma_Estimate","True_TF_Prop")
  return(output)
}


#4.
##Plotting function for holding number of trees constant.
##Takes output from bart_loop()
bart_plot_const_tree=function(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent){
  ##Initial Plot Set-up
    par(mgp=c(1.8,.5,0))
  col.vec=c("orange","green","blue","red")
  plot(1,type="n",xlim=c(min(ntree.vec),max(ntree.vec)+4),ylim=c(min(factor.vec),max(factor.vec)+1),
       xlab="Number of Trees",ylab="Noise/Signal Ratio",
       main=paste("Break Levels for BART Models\n",num.obs,"Observations"))
  
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
bart_plot_const_tf=function(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent){
  ##Initial Plot Setup
    par(mgp=c(1.8,.5,0))
  col.vec=c("orange","green","blue","red")
  plot(1,type="n",ylim=c(min(factor.vec),max(factor.vec)+1),xlim=c(min(tfsize.vec),max(tfsize.vec)+8),
       xlab="Number of TFs",ylab="Noise/Signal Ratio",
       main=paste("Break Levels for BART Models\n",num.obs,"Observations"),
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
bart_plot_both=function(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent){
  par(mfrow=c(1,2))
  bart_plot_const_tree(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent)
    bart_plot_const_tf(output,ntree.vec,factor.vec,tfsize.vec,num.obs, percent)
}


#7.
##Plotting-4 plots per forest size. each based on number of tfs. 
##Each plot contains proportions plotted by N/S ratio. Then 1 for each of 4 types
bart_loop_plot=function(which.treesize,output,tf.size.vec){
  par(mfrow=c(2,2))
  par(mgp=c(1.8,.5,0))
  for(i in 1:length(tf.size.vec)){
    temp=subset(output,subset=Num_Trees==which.treesize & Num_TFs==tf.size.vec[i] )
  
  
    plot(temp$"Noise/Signal",temp$"1st",type="o",lwd=2,lty=2,col="green",
        main=paste("Split Selection Proportion\n Num_TFs=",tf.size.vec[i],",Num_Trees=",which.treesize)
        ,xlab="Noise-To-Signal Ratio",ylab="Proportion")
    lines(temp$"Noise/Signal",temp$"2nd",type="o",lwd=3,lty=2,col="blue")
    lines(temp$"Noise/Signal",temp$"3rd",type="o",lwd=2,lty=2,col="red")
    lines(temp$"Noise/Signal",temp[,14],type="o",lwd=2,lty=2,col="black")
#legend("topright",c("True TF","1st Selected","2nd Selected","3rd Selected"),
 #        lwd=2,lty=2,col=c("black","green","blue","red"),cex=.85)
    print(i)
    }
  }

##Aggregating Functions for Plotting
##8.
tf_axis_aggregate=function(output,tfsize.vec,ntree.vec,percent,fn_name){
  output.vec=numeric(length(tfsize.vec))
  fun=match.fun(fn_name)
  for(i in 1:length(tfsize.vec)){
    temp.vec.for.calc=numeric(length(ntree.vec))
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
      temp.vec.for.calc[j]=temp$"Noise/Signal"[noise_count]   
    }
    output.vec[i]=median(temp.vec.for.calc)
    temp.vec.for.calc=numeric(length(ntree.vec))
  }
  return(output.vec)
}


##9
tree_axis_aggregate=function(output,tfsize.vec,ntree.vec,percent,fn_name){
  output.vec=numeric(length(ntree.vec))
  fun=match.fun(fn_name)
  for(i in 1:length(ntree.vec)){
    temp.vec.for.calc=numeric(length(tfsize.vec))
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
      temp.vec.for.calc[j]=temp$"Noise/Signal"[noise_count]   
    }
    output.vec[i]=median(temp.vec.for.calc)
    temp.vec.for.calc=numeric(length(tfsize.vec))
  }
  return(output.vec)
}
