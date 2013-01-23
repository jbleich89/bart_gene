#!/bin/bash

#$ -m e -M kapelner@wharton.upenn.edu
#$ -j y
#$ -N tf_discovery_for_gene_expression
#$ -t 2

let START=$END
let END=$SGE_TASK_ID*1

echo "starting R task # $SGE_TASK_ID START: $START END: $END"

R --no-save '--args GENE_START=$START GENE_END=$END' < analysis_with_tuning/simulation.R