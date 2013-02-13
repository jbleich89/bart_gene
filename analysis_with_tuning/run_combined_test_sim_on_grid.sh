#!/bin/bash

#$ -j y
#$ -N tf_discovery_for_gene_expression_combined_tests
#$ -t 1-5000
#$ -q intel

echo "starting combined test for gene # $SGE_TASK_ID"
export _JAVA_OPTIONS="-Xms128m -Xmx5000m"

R --no-save --args iter_num=$SGE_TASK_ID < analysis_with_tuning/combined_tests.R