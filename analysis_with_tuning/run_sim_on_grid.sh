#!/bin/bash

#$ -m e -M kapelner@wharton.upenn.edu
#$ -j y

R --no-save '--args GENE_START=1 GENE_END=2' < analysis_with_tuning/simulation.R