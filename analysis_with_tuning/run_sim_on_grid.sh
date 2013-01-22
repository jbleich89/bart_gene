#!/bin/bash

# This will execute R in our batch mode
R CMD BATCH --no-save --no-restore '--args GENE_START=1 GENE_END=10' analysis_with_tuning/simulation.R simulation.out &

#qsub -m e -M username@wharton.upenn.edu your-job-script.sh

#echo 'hostname; echo $SGE_TASK_ID; sleep 20' | qsub -N ArrayTest1 -t 1-5