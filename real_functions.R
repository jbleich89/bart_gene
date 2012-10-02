#Functions for real data


##Functions###########
##Rounding function for setting prior
roundUpNice <- function(x, nice=c(1,2,3,4,5,6,7,8,9,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

changePrior=function(location){
  temp=location*10
  ifelse(temp>1,roundUpNice(temp),1)
}

ChangePriorVec=function(locVec){
  sapply(locVec,changePrior)
}
