##SVD Script 

N = 250
p = 40
NUM_DATA = 100
NUM_BART = 50
SIGMA = 1
BURN_SIZE = 2000
POST_SIZE = 5000
NTREE = 10

probs_mat_svd = matrix(NA, nrow = NUM_DATA*NUM_BART, ncol = p + 2)
corrs_mat_svd = matrix(NA, nrow = NUM_DATA, ncol = p)
counter = 1

for(i in 1 : NUM_DATA){
  x.train = data.frame(matrix(rnorm(N * p,0,1),nrow = N, ncol = p))
  y = rnorm(N, 0, SIGMA)
  decomp = svd(cbind(y,x.train))
  u = decomp$u
  dim(cor(u))
  uy = u[, 1]
  ux = data.frame(u[, 2 : (p + 1)])
  ##get correlations
  corrs_mat_svd[i, ] = sapply(1 : p, function(s) abs(cor(y, x.train[,s])))
  for(j in 1 : NUM_BART){
    ## do bart stuff
    bart_noise = build_bart_machine(X = ux, y = uy, num_trees = NTREE, num_burn_in = BURN_SIZE, num_iterations_after_burn_in = POST_SIZE, verbose = F, run_in_sample = F)
    avg_props = get_var_props_over_chain(bart_noise)
    destroy_bart_machine(bart_noise)
    
    probs_mat_svd[counter, ] = c(i, j , avg_props)
    
    if(counter %% 25 == 0) print(counter)
    counter = counter + 1
  }
}

save(probs_mat_svd, file = "null_probs_matrix_svd.Rdata")
save(corrs_mat_svd, file = "null_corrs_matrix_svd.Rdata")