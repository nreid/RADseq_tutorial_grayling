#!/bin/bash 
#SBATCH --job-name=freebayes
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 30
#SBATCH --mem=90G
#SBATCH --qos=general
#SBATCH --partition=general

#load software
module load htslib/1.10.2
module load freebayes/1.3.1

#input, output files, directories
INDIR=../results/aligned/
OUTDIR=../results/freebayes
mkdir -p $OUTDIR

#reference genome
GEN=../genome/GCA_004348285.1_ASM434828v1_genomic.fna

# make a list of bam files
	# use popmap_total.txt to construct the list
POPMAP=../meta/popmap_total.txt

BAMLIST=$OUTDIR/bam.list
cat $POPMAP | cut -f 1 | sed 's/$/.bam/' | sed  "s,^,$INDIR," >$BAMLIST

./freebayes-parallel \
	<(python fasta_generate_regions.py ${GEN}.fai 5000000) 30 \
	-f ${GEN} \
	--bam-list $OUTDIR/bam.list \
	-m 30 \
	-q 20 \
	--min-coverage 1500 \
	--skip-coverage 30000 | \
bgzip -c >$OUTDIR/fb_parallel.vcf.gz

# index the vcf file
tabix -p vcf $OUTDIR/fb_parallel.vcf.gz

