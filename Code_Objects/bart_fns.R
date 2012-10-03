##BART Functions s
##calculate variable sums
sum_calc=function(bart.list){apply(bart.list[[3]],2,sum)}

##calculate variable proportions 
prop_calc=function(bart.list){
  tot=sum(apply(bart.list[[3]],2,sum))
  props=apply(bart.list[[3]],2,sum)/tot
  round(props,3)
  }

##Plotting function for variable selection
##must load in order 50,20,10 
plot_vars=function(full.bart,num.true.tf,num.tf,title){
  prop50=prop_calc(full.bart[[1]])
  prop20=prop_calc(full.bart[[2]])  
  prop10=prop_calc(full.bart[[3]])
  y.min=min(min(prop50),min(prop20),min(prop10))
  y.max=max(max(prop50),max(prop20),max(prop10))
  par(mgp=c(1.8,.5,0))
  plot(1:num.tf,prop50,type="o",
       main=title,
       xlab="TF Number",
       ylab="Inclusion Proportion",
       ylim=c(y.min,y.max),
       lty=2,
       pch=21,
       col="red",
       xaxt='n')
  box()
  axis(side=1,at=2*(1:num.tf))
  lines(1:num.tf,prop20,type="o",lty=3,lwd=1,col="blue",pch=22)
  lines(1:num.tf,prop10,type="o",lty=4,col="forestgreen",pch=23)
  abline(v=num.true.tf,lty=1,lwd=1,col="black")
  legend("topright",c("50 Trees","20 Trees","10 Trees"),
         pch=21:23,lty=2:4,col=c("red","blue","forestgreen"))
}

plot_vars=function(full.bart,num.true.tf,num.tf,title, thin){
  prop50=prop_calc(full.bart[[1]])
  prop20=prop_calc(full.bart[[2]])  
  prop10=prop_calc(full.bart[[3]])
  y.min=min(min(prop50),min(prop20),min(prop10))
  y.max=max(max(prop50),max(prop20),max(prop10))
  par(mgp=c(1.8,.5,0))
  plot(1:num.tf,prop50,type="o",
       main=title,
       xlab="TF Number",
       ylab="Inclusion Proportion",
       ylim=c(y.min,y.max),
       lty=2,
       pch=21,
       col="red",
       xaxt='n')
  box()
  axis(side=1,at=2*(1:num.tf))
  lines(1:num.tf,prop20,type="o",lty=3,lwd=1,col="blue",pch=22)
  lines(1:num.tf,prop10,type="o",lty=4,col="forestgreen",pch=23)
  abline(v=num.true.tf,lty=1,lwd=1,col="black")
  legend("topright",c("50 Trees","20 Trees","10 Trees"),
         pch=21:23,lty=2:4,col=c("red","blue","forestgreen"))
}

plot_vars_thin=function(full.bart,num.true.tf,num.tf,title,thin){
  
  prop_calc=function(bart.obj){
  tot=sum(apply(bart.obj[[3]][thin*(1:5000/thin),],2,sum))
  props=apply(bart.obj[[3]][thin*(1:5000/thin),],2,sum)/tot
  }
  
  prop50=prop_calc(full.bart[[1]])
  prop20=prop_calc(full.bart[[2]])  
  prop10=prop_calc(full.bart[[3]])
  y.min=min(min(prop50),min(prop20),min(prop10))
  y.max=max(max(prop50),max(prop20),max(prop10))
  par(mgp=c(1.8,.5,0))
  plot(1:num.tf,prop50,type="o",
       main=title,
       xlab="TF Number",
       ylab="Inclusion Proportion",
       ylim=c(y.min,y.max),
       lty=2,
       pch=21,
       col="red",
       xaxt='n')
  box()
  axis(side=1,at=2*(1:num.tf))
  lines(1:num.tf,prop20,type="o",lty=3,lwd=1,col="blue",pch=22)
  lines(1:num.tf,prop10,type="o",lty=4,col="forestgreen",pch=23)
  abline(v=num.true.tf,lty=1,lwd=1,col="black")
  legend("topright",c("50 Trees","20 Trees","10 Trees"),
         pch=21:23,lty=2:4,col=c("red","blue","forestgreen"))
}


##plotting function for burn in
plot_sig=function(bart.list,n_burn,subtitle){
  sig=bart.list[[2]]
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

num_correct=function(bart.list,num.tf.true){
  
  prop_calc=function(bart.obj){
  tot=sum(apply(bart.obj[[3]],2,sum))
  props=apply(bart.obj[[3]],2,sum)/tot
  }
  props=prop_calc(bart.list)
  sorted.props=sort(x=props,decreasing=T)
  correct=intersect(props[1:num.tf.true],sorted.props[1:num.tf.true])
  return(length(correct))
}


num_correct_thin=function(bart.list,num.tf.true,thin){
  
  prop_calc=function(bart.obj){
  tot=sum(apply(bart.obj[[3]][thin*(1:5000/thin),],2,sum))
  props=apply(bart.obj[[3]][thin*(1:5000/thin),],2,sum)/tot
  }
  props=prop_calc(bart.list)
  sorted.props=sort(x=props,decreasing=T)
  correct=intersect(props[1:num.tf.true],sorted.props[1:num.tf.true])
  return(length(correct))
}