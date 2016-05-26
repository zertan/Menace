#!/usr/bin/env bash
# Note: these paths and variable names are exact and must
# lead to existing directories and files organized according
# to the directory view @ https://github.com/zertan/PTR-Pipeline.

export DATA_PATH=""					# path to metagenomic cohort (top directory)
export REF_PATH=""					# the path to the bowtie2 index
export REFNAME=""					# the name of the bowtie2 index $REFNAME.*.bt2l
export OUT_PATH=""					# ouput path of cluster data	
export DORIC_PATH=""				# 
export EMAIL=""						# email to be used by slurm and fetchSeq.py
export SUBMIT_PATH="$(pwd)"			# path to PTR-Build where you execute slurm (sbatch)
export CLUSTER=""					# the name of the cluster
export CPUCORES=""					# the number of cores per compute node

