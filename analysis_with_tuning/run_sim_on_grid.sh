#!/bin/bash

#$ -m e -M kapelner@wharton.upenn.edu
#$ -j y
# This will execute R in our batch mode

R --no-save '--args GENE_START=1 GENE_END=10' analysis_with_tuning/simulation.R