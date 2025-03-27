#!/bin/bash

#SBATCH --job-name=med_cc_${1}_${2}
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --time=3-23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=epwkb
#SBATCH --export=ALL
#SBATCH --account=SSCM026803
### SBATCH --partition=short
#SBATCH --array=1

# Load in Stata from the module list
module add apps/stata/17

# CD
cd /user/home/epwkb/scratch/inflam/complete_case/${1}/${2}

# Application name
application1="stata do /user/home/epwkb/scratch/inflam/complete_case/stata_imps_array.do"

# Array job
start_time=`date +%s`

echo ${SLURM_ARRAY_TASK_ID}
echo ${1}
echo ${2}

# Application
/usr/bin/time -v $application1 ${SLURM_ARRAY_TASK_ID} ${1} ${2}

end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.

###
