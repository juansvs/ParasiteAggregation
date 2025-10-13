#!/bin/bash

#SBATCH --job-name=FoxSim
#SBATCH --output=FoxSim_%J.out
#SBATCH --error=FoxSim_%J.err
#SBATCH --time=24:00:00
#SBATCH --partition=k2-medpri
#SBATCH --ntasks=10
#SBATCH --nodes=4
#SBATCH --mem-per-cpu=2G

#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.vargas@qub.ac.uk

module load apps/R/4.4.1/gcc-14.1.0+openblas-0.3.27


cd ~/ParasiteAggregation
Rscript run_FoxIBM_parallel.R
