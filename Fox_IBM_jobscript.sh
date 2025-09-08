#!/bin/bash

#SBATCH --job-name=FoxSim
#SBATCH --output=FoxSim-output-file
#SBATCH --error=FoxSim-error-file
#SBATCH --time=24:00:00
#SBATCH --partition=k2-medpri
#SBATCH --ntasks=50
#SBATCH --nodes=1-10

#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.vargas@qub.ac.uk

module load apps/R/4.4.1/gcc-14.1.0+openblas-0.3.27


cd ~/ParasiteAggregation
Rscript run_FoxIBM_parallel.R