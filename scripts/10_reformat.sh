#!/bin/bash
#SBATCH --job-name=reformat
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bcftools/1.9
module load plink/2.00a2.3LM


# input/output directories
INDIR=../results/filtered_vcfs
OUTDIR=../results/variants_reformatted
mkdir -p $OUTDIR

# generate input file(s) for admixture
plink2 --vcf $INDIR/fb_final.vcf.gz --allow-extra-chr --out $OUTDIR/fb

# create a reformatted file with genotypes only
bcftools query -H -f '%CHROM\t%POS\t[\t%GT]\n' $INDIR/fb_final.vcf.gz >$OUTDIR/fb_genos.tsv

