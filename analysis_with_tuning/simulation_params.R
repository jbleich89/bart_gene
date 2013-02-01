#simulation params

if (NOT_ON_GRID){
	cs = as.character(c(0))	

} else {
	cs = as.character(c(
		0, #no prior
		0.5, #weak prior
		1,
		2,
		4, #medium strength prior
		10000 #prior dominates
	))	
}


ms = as.character(c(
#	10,
	20
))

alphas = as.character(c(
#	0.01,
	0.05
#	0.10
))

if (NOT_ON_GRID){
	NUM_PERMUTE_SAMPLES = 3
	NUM_REP_FOR_TRAIN = 2
} else {
	NUM_PERMUTE_SAMPLES = 100
	NUM_REP_FOR_TRAIN = 10
}

NUM_BURN_IN = 2000
NUM_ITER_AFTER = 2000

NUM_CORES = 1
ALPHA = 0.05
METHODS = c("important_tfs_at_alpha_pointwise", "important_tfs_at_alpha_simul_max", "important_tfs_at_alpha_simul_se")
NUM_TREES_FOR_EVAL = 200
GENE_NUM = 1 #this will be set by passing an arg into R below:

source("helper_functions.R")

if (NOT_ON_GRID){
	setwd(paste("C:/Users/", LAST_NAME, "/Desktop/Dropbox/BART_gene", sep = ""))
} else {
	setwd("../../shane_data")
}

priors = read.table("CHIP.priorprobs.39.txt", header = TRUE)
gene.exp = read.table("expression.genes.txt", header = TRUE)
tf.exp = read.table("expression.tfs.39.txt", header = TRUE)

priorMat = as.matrix(priors[, 3 : ncol(priors)])
rownames(priorMat) = priors[, 1] 
colnames(priorMat) = colnames(priors)[3 : ncol(priors)]
prior_mat_adj = prior_adj(priorMat)
result = setup() ##5 objects
gene_train = result[["gene.train"]]
gene_cv = result[["gene.cv"]]
gene_test = result[["gene.test"]]##gene.train and gene.test rows are obs and cols are genes
tf_train = result[["tf.train"]]
tf_cv = result[["tf.cv"]]
tf_test = result[["tf.test"]]##TF train and TF test: rows are obs and cols are TFs
gene_names = as.character(gene.exp[, 2]) ##gene names
