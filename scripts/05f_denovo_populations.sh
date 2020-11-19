#!/bin/bash 
#SBATCH --job-name=populations
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load stacks/2.53
module load GATK/4.1.8.1

# input, output files, directories

INDIR=../results/stacks/denovo

POPMAP=../meta/popmap_total.txt

populations \
-P $INDIR \
-M $POPMAP \
-p 1 \
-R 0.50 \
--hwe \
--genepop \
--vcf \
--fasta-samples \
--fasta-loci \
--treemix \
--structure \
-t 8
