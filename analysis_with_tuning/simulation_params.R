#simulation params

cs = c(
	0, #no prior
	0.5, #weak prior
	1,
	2,
	4, #medium strength prior
	10000 #prior dominates
)

ms = c(
#	10,
	20
)

NUM_PERMUTE_SAMPLES = 2
NUM_BURN_IN = 2000
NUM_ITER_AFTER = 2000
NUM_REP_FOR_TRAIN = 1
NUM_CORES = 1
ALPHA = 0.05
METHODS = c("important_tfs_at_alpha_pointwise", "important_tfs_at_alpha_simul_max", "important_tfs_at_alpha_simul_se")
NUM_TREES_FOR_EVAL = 200
GENE_NUM = 1 #this will be set by passing an arg into R below:
