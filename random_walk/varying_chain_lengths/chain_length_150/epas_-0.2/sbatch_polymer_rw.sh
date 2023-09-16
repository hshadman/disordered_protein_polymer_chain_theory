#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --job-name=15_-0.2
#SBATCH --partition=computeq
#SBATCH --nodes=1
#SBATCH --time=28-00:00:00
#SBATCH --mem=20G
##SBATCH --error=job.%J.err
##SBATCH --output=job.%J.out

cd $SLURM_SUBMIT_DIR
echo "$SLURM_JOB_ID"
#run my job below 
#check directory I am in

pwd
echo "Started at time" 
echo date

# make sure you have initial.dat and 44rep.input prepared in the directory that you are issuing this sbatch
#make sure you point to the right location of the executable relative to the directory you are submitting 

../../../RWchain_v0

echo "finished RWchain_v0  at"
date






