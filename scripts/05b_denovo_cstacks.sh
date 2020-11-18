#!/bin/bash 
#SBATCH --job-name=cstacks
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

# load software
module load stacks/2.53
module load R/3.6.3

# input, output files, directories
INDIR=../results/stacks/denovo

# create popmap file
	# this R script creates the popmap files required for the next 
	# steps of stacks. we're using an R script so we can subset 
	# samples for cstacks, which would take too long with the total 
	# dataset. the next line of code executes an R script that 
	# does this. it's commented out for the tutorial because it 
	# requires the tidyverse packages to be installed, which users 
	# may not have, but is included to completely document steps required
	# 
	# 
# Rscript make_popmaps.R

POPMAP=../meta/popmap_cstacks.txt

awk '$13 ~ /LandscapeGenomics/' $METADATA | cut -f 3,4 >$POPMAP

cstacks \
-P $INDIR \
-M $POPMAP \
-n 2
-p 10 