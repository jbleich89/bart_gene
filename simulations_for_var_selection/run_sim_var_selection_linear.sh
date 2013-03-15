#!/bin/bash

#$ -j y
#$ -N bart_for_variable_selection_linear_sims
#$ -t 1-125
#$ -q intel

echo "starting R for gene # $SGE_TASK_ID"
export _JAVA_OPTIONS="-Xms128m -Xmx5300m"

R --no-save --args iter_num=$SGE_TASK_ID < simulations_for_var_selection/sims_for_bart_var_selection_lin_mod.R