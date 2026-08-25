#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p veryhimem
#SBATCH -c 2
#SBATCH --mem=66000M
#SBATCH -t 0-20:00 # Runtime in D-HH:MM
#SBATCH -J run_Rscript_all

work_path=$SLURM_SUBMIT_DIR
cd $work_path

module load R/4.1.0

Rscript  DiffBind_analysis_ATLL_all_comp.R
#Rscript DiffBind_analysis_HBL1_all_comp.R

