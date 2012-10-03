########Aggregate Plotting Functions################

####TF Axis
bart_agg_plot_tf_axis=function(tf.const.mat,tfsize.vec,factor.vec){
  par(mgp=c(1.8,.5,0),mar=c(3,3,3,2))
  col.vec=c("orange","green","blue","red","black")
  plot(1,type="n",xlim=c(min(tfsize.vec),max(tfsize.vec)+6),ylim=c(min(factor.vec),max(factor.vec)+1),
       xlab="Number of TFs",ylab="Noise/Signal Ratio",
       main=paste("Break Levels for BART Models\n"),xaxt="n")
  axis(side=1,at=tfsize.vec)
  for(i in 1:nrow(tf.const.mat)){
    lines(tfsize.vec,tf.const.mat[i,2:ncol(tf.const.mat)],type="o",col=col.vec[i],pch=16,lwd=2)
  }
  leg.names=paste(tf.const.mat[,1],"Obs.")
  legend("topright",legend=leg.names,lwd="2",pch=16,col=col.vec,cex=.7)
}

##########3
tf_axis_aggregate=function(output,tfsize.vec,ntree.vec,percent,fn_name){
  output.vec=numeric(length(tfsize.vec))
  fun=match.fun(fn_name)
  print(fun)
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
    print(temp.vec.for.calc)
    output.vec[i]=fun(temp.vec.for.calc)
    temp.vec.for.calc=numeric(length(ntree.vec))
  }
  return(output.vec)
}



#######Tree Axis###########3
bart_agg_plot_tree_axis=function(tree.const.mat,ntree.vec,factor.vec){
  par(mgp=c(1.8,.5,0),mar=c(3,3,3,2))
  col.vec=c("orange","green","blue","red","black")
  plot(1,type="n",xlim=c(min(ntree.vec),max(ntree.vec)+4),ylim=c(min(factor.vec),max(factor.vec)+1),
       xlab="Number of Trees",ylab="Noise/Signal Ratio",
       main=paste("Break Levels for BART Models\n"))
  for(i in 1:nrow(tree.const.mat)){
    lines(ntree.vec,tree.const.mat[i,2:ncol(tree.const.mat)],type="o",col=col.vec[i],pch=16,lwd=2)
  }
  leg.names=paste(tree.const.mat[,1],"Obs.")
  legend("topright",legend=leg.names,lwd="2",pch=16,col=col.vec,cex=.7)
}

############
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
    output.vec[i]=fun(temp.vec.for.calc)
    temp.vec.for.calc=numeric(length(tfsize.vec))
  }
  return(output.vec)
}
