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

# input, output files, directories
INDIR=../results/stacks/denovo

# create popmap file
# we'll only use samples labeled as part of the landscape genomics subsection of the project
METADATA=../meta/Urban_ddRAD_FishIDs_Bioinformatics_2018.tsv
POPMAP=../meta/popmap.txt

awk '$13 ~ /LandscapeGenomics/' $METADATA | cut -f 3,4 >$POPMAP

cstacks \
-P $INDIR \
-M $POPMAP \
-n 2
-p 10 