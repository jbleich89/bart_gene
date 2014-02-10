##Dynatree variable selection
var_sel_dynaTree = function(X, y,  n_particles = 1000){ ##package default
  selected_col_names = colnames(X)
  ##do first round and select some variables
  dt = dynaTree(X = as.matrix(X), y = y, model = "constant", N = n_particles, verb=0)
  dt = relevance(dt, approx = T)
  
  ##select variables using their metric -- postive value
  selected_col_names = selected_col_names[which(apply(dt$relevance, 2 ,median) > 0)]
  
  ##now loop through backward selection until the Bayes Factor is smaller
  bf = -1
  while(bf < 0){
    dt_prop = dynaTree(X = as.matrix(X[ ,selected_col_names]), y = y, model = "constant", N = n_particles, verb=0)
    bf = mean(getBF(dt, dt_prop)) ##using mean bayes factor seems okay to me
    #print(bf)
    if(bf < 0){ ##this means log(M1) - log(M2) < 0 so the posterior for M1 is less than M2 so choose M2. 
      dt = dt_prop
      dt = relevance(dt, approx = T)
      head(dt$relevance)
      selected_col_names_new = selected_col_names[which(apply(dt$relevance, 2 ,median) > 0)]
      if(length(selected_col_names) == length(selected_col_names_new)){
        if(all(selected_col_names_new == selected_col_names)) return(selected_col_names) ##break out if you get the same thing.
        selected_col_names = selected_col_names_new
      }
    }
  }
  selected_col_names
}
