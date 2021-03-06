##Read partial results

setwd("C:/Users/jbleich/Dropbox/BART_gene/sim_results_new/linear/")

#p200 p0 20 sigsq 5
results_p200_p20_sigsq5 = array(dim= c(20, 50, 3))
dimnames(results_p200_p20_sigsq5) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p200_p0_20_sigsq_5_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p200_p0_20_sigsq_5_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p200_p20_sigsq5[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p200_p20_sigsq5, file = "results_p200_p20_sigsq5.Rdata")



#p200 p0 20 sigsq 20
results_p200_p20_sigsq20 = array(dim= c(20, 50, 3))
dimnames(results_p200_p20_sigsq20) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p200_p0_20_sigsq_20_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p200_p0_20_sigsq_20_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p200_p20_sigsq20[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p200_p20_sigsq20, file = "results_p200_p20_sigsq20.Rdata")


#p200 p0 2 sigsq 5
results_p200_p2_sigsq5 = array(dim= c(20, 50, 3))
dimnames(results_p200_p2_sigsq5) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p200_p0_2_sigsq_5_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p200_p0_2_sigsq_5_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p200_p2_sigsq5[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p200_p2_sigsq5, file = "results_p200_p2_sigsq5.Rdata")



#p200 p0 2 sigsq 20
results_p200_p2_sigsq20 = array(dim= c(20, 50, 3))
dimnames(results_p200_p2_sigsq20) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p200_p0_2_sigsq_20_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p200_p0_2_sigsq_20_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p200_p2_sigsq20[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p200_p2_sigsq20, file = "results_p200_p2_sigsq20.Rdata")



#p500 p0 25 sigsq 1
results_p500_p25_sigsq1 = array(dim= c(20, 50, 3))
dimnames(results_p500_p25_sigsq1) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p500_p0_25_sigsq_1_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p500_p0_25_sigsq_1_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p500_p25_sigsq1[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p500_p25_sigsq1, file = "results_p500_p25_sigsq1.Rdata")



#p500 p0 25 sigsq 5
results_p500_p25_sigsq5 = array(dim= c(20, 50, 3))
dimnames(results_p500_p25_sigsq5) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p500_p0_25_sigsq_5_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p500_p0_25_sigsq_5_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p500_p25_sigsq5[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p500_p25_sigsq5, file = "results_p500_p25_sigsq5.Rdata")


#p500 p0 50 sigsq 1
results_p500_p50_sigsq1 = array(dim= c(20, 50, 3))
dimnames(results_p500_p50_sigsq1) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p500_p0_50_sigsq_1_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p500_p0_50_sigsq_1_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p500_p50_sigsq1[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p500_p50_sigsq1, file = "results_p500_p50_sigsq1.Rdata")



#p500 p0 50 sigsq 5
results_p500_p50_sigsq5 = array(dim= c(20, 50, 3))
dimnames(results_p500_p50_sigsq5) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_linear_p500_p0_50_sigsq_5_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_linear_p500_p0_50_sigsq_5_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_p500_p50_sigsq5[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_p500_p50_sigsq5, file = "results_p500_p50_sigsq5")

##

apply(results_p500_p50_sigsq1[,,3], 1, mean)  ##gets average F1
apply(results_p200_p20_sigsq5[,,3], 1, sd)/5  ##gets average F1

apply(results_fried_unif_p500_sigsq1[,,2], 1, mean)


##Friedman Uniform
setwd("C:/Users/jbleich/Dropbox/BART_gene/sim_results_new/friedman_uniform/")

#p500 sigsq 1
results_fried_unif_p500_sigsq1 = array(dim= c(20, 50, 3))
dimnames(results_fried_unif_p500_sigsq1) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                           "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                           "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                           "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_friedman_uniform_p500_sigsq_1_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_friedman_uniform_p500_sigsq_1_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_fried_unif_p500_sigsq1[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_fried_unif_p500_sigsq1, file = "results_fried_unif_p500_sigsq1.Rdata")

#p500 sigsq 25
results_fried_unif_p500_sigsq25 = array(dim= c(20, 50, 3))
dimnames(results_fried_unif_p500_sigsq25) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                                  "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                                  "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                                  "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_friedman_uniform_p500_sigsq_25_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_friedman_uniform_p500_sigsq_25_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_fried_unif_p500_sigsq25[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_fried_unif_p500_sigsq25, file = "results_fried_unif_p500_sigsq25.Rdata")



#p500 sigsq 100
results_fried_unif_p500_sigsq100 = array(dim= c(20, 50, 3))
dimnames(results_fried_unif_p500_sigsq100) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                                  "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                                  "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                                  "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_friedman_uniform_p500_sigsq_100_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_friedman_uniform_p500_sigsq_100_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_fried_unif_p500_sigsq100[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_fried_unif_p500_sigsq100, file = "results_fried_unif_p500_sigsq100.Rdata")


setwd("C:/Users/jbleich/Dropbox/BART_gene/sim_results_new/friedman_uniform/")

#p1000 sigsq 1
results_fried_unif_p1000_sigsq1 = array(dim= c(20, 50, 3))
dimnames(results_fried_unif_p1000_sigsq1) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                                  "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                                  "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                                  "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_friedman_uniform_p1000_sigsq_1_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_friedman_uniform_p1000_sigsq_1_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_fried_unif_p1000_sigsq1[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_fried_unif_p1000_sigsq1, file = "results_fried_unif_p1000_sigsq1.Rdata")

#p1000 sigsq 25
results_fried_unif_p1000_sigsq25 = array(dim= c(20, 50, 3))
dimnames(results_fried_unif_p1000_sigsq25) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                                   "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                                   "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                                   "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_friedman_uniform_p1000_sigsq_25_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_friedman_uniform_p1000_sigsq_25_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_fried_unif_p1000_sigsq25[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_fried_unif_p1000_sigsq25, file = "results_fried_unif_p1000_sigsq25.Rdata")



#p1000 sigsq 100
results_fried_unif_p1000_sigsq100 = array(dim= c(20, 50, 3))
dimnames(results_fried_unif_p1000_sigsq100) = list(c("BART_CV", "BART_pointwise", "BART_simul_max", "BART_simul_se", "BART_CV_good_prior", "BART_pointwise_good_prior", 
                                                    "BART_simul_max_good_prior", "BART_simul_se_good_prior", "BART_CV_bad_prior", "BART_pointwise_bad_prior", 
                                                    "BART_simul_max_bad_prior",  "BART_simul_se_bad_prior", "stepwise_backward", "stepwise_forward", 
                                                    "lasso", "RF_CV", "RF_point", "RF_simul", "dynaTree","spikeSlab") , NULL, c("Precision", "Recall", "F1"))
for(nr in 1:50){
  temp = read.csv(paste("partial_var_sel_sim_friedman_uniform_p1000_sigsq_100_nr_", nr, "_revised.csv", sep = ""))
  temp2 = read.csv(paste("partial_var_sel_sim_friedman_uniform_p1000_sigsq_100_nr_", nr, "_bart_gmax.csv", sep = ""))
  temp[c(3,7, 11),2:3] = temp2[c(3,7, 11),2:3]
  f1s = 2 * temp[,2] * temp[,3]/(temp[,2] + temp[, 3])
  f1s = ifelse(is.nan(f1s), 0, f1s)
  results_fried_unif_p1000_sigsq100[,nr,] = cbind(temp[,2], temp[,3], f1s)
}
save(results_fried_unif_p1000_sigsq100, file = "results_fried_unif_p1000_sigsq100.Rdata")

