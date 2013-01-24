#!/bin/bash

#$ -m e -M kapelner@wharton.upenn.edu
#$ -j y
#$ -N tf_discovery_for_gene_expression
#$ -t 1-20

echo "starting R for gene # $SGE_TASK_ID"
export _JAVA_OPTIONS="-Xms128m -Xmx3g"

R --no-save --args GENE_NUM=$SGE_TASK_ID < analysis_with_tuning/simulation.R