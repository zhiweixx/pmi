#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 100G
#SBATCH -t 36:00:00
#SBATCH -p medium
#SBATCH -o err/job_%j.out 

# Execute commands
module load gcc/6.2.0 R/4.0.1
R CMD BATCH "--args $arg1 $arg2 $arg3 $arg4" SIM_ON_THMS4.R err/moon.out


