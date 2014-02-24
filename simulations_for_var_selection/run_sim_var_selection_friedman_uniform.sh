#!/bin/bash

#$ -j y
#$ -N bart_for_variable_selection_friedman_sims_uniform
#$ -t 1-6
#$ -q intel

echo "starting R for task # $SGE_TASK_ID"
export _JAVA_OPTIONS="-Xms128m -Xmx5300m"

R --no-save --args iter_num=$SGE_TASK_ID < simulations_for_var_selection/sims_for_bart_var_selection_friedman_mod_uniform.R
