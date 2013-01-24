#!/bin/bash

#$ -j y
#$ -N tf_discovery_for_gene_expression
#$ -t 1-100
#$ -q bigram

echo "starting R for gene # $SGE_TASK_ID"
export _JAVA_OPTIONS="-Xms128m -Xmx7000m"

R --no-save --args GENE_NUM=$SGE_TASK_ID < analysis_with_tuning/simulation.R