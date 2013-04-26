tryCatch(library(energy), error = function(e){install.packages("energy")}, finally = library(energy))
tryCatch(library(xtable), error = function(e){install.packages("xtable")}, finally = library(xtable))

LAST_NAME = "kapelner"
NOT_ON_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

#######boxplots
if (.Platform$OS.type == "windows"){
	setwd("C:/Users/Kapelner/workspace/bart_gene/analysis_with_tuning")
}

source("simulation_params.R")

setwd("C:/Users/Kapelner/Desktop/Dropbox/BART_gene")
load("all_results.RData")
load("all_validations.RData")
load("all_rmse_results.RData")
load("all_num_var_results.RData")


MAX_GENE_NUM = 6026