#!/bin/bash

#SBATCH --job-name=imputations_31102023
#SBATCH -o imputations_31102023-output
#SBATCH -e imputations_31102023-error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=2-15:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=epwkb
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --account=SSCM026803

# Modules
#module add languages/R/4.3.3

#interactive: srun --ntasks=1  --cpus-per-task=4 --time=24:00:00 --partition=test --pty bash -I # for interactive session on node

# CD
cd /user/work/epwkb/aces/scripts/02-ACEs/

#! Application name
application1="Rscript 03_Script_3_HPC.R"

#! Run options
#options=""

#####################################
#####################################

echo Imputation test run: SLURM edition

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job ID is $SLURM_JOBID
echo This jobs runs on the following machines:
echo $SLURM_JOB_NODELIST

#! Run the threaded exe
date
#srun /usr/bin/time -v $application1 $options
/usr/bin/time -v $application1 $options
date
