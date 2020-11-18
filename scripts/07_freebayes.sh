#!/bin/bash 
#SBATCH --job-name=freebayes
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 7
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


hostname
date

# load required software

module load htslib/1.10.2
module load freebayes/1.3.1

# input, output directories
INDIR=../results/aligned/

OUTDIR=../results/freebayes
mkdir -p $OUTDIR 

# make a list of bam files
	# use popmap_total.txt to construct the list
POPMAP=../meta/popmap_total.txt

BAMLIST=$OUTDIR/bam.list
cat $POPMAP | cut -f 1 | sed 's/$/.bam/' | sed  "s,^,$INDIR," >$BAMLIST

# set a variable for the reference genome location
GEN=../genome/GCA_004348285.1_ASM434828v1_genomic.fna

# run freebayes
	# skip regions with excess or deficient coverage
freebayes \
-f $GEN \
--bam-list $OUTDIR/bam.list \
-m 30 \
-q 20 \
--min-coverage 1500 \
--skip-coverage 30000 | \
bgzip -c >$OUTDIR/fb.vcf.gz

# index the vcf file
tabix -p vcf $OUTDIR/fb.vcf.gz

date