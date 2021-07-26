#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2g
#SBATCH -t 5-00:00:00

Rscript R_code/R_simulations/Sim_1.R

