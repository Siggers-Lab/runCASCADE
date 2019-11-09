#!/bin/bash

#$ -N CASCADE_a4_run1
#$ -cwd
#$ -j y
#$ -V
#$ -m be
#$ -M ykoga07@bu.edu
#$ -l h_rt=12:00:00

module load R/3.5.1
Rscript CASCADE_analysis.R a4_run1 ENH_CASCADE_040119 ENH_TILE_001_FULL_ANNOT.bed /projectnb/siggers/data/yusuke_project/PBM_matrices/ a4_run1_exp.txt ENH_TILE_001_TF_sites.txt TRUE
