#!/bin/bash
#SBATCH --job-name=filter_vcfs
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
module load htslib/1.9
module load vcftools/0.1.16
module load vcflib/1.0.0-rc1

OUTDIR=../results/filtered_vcfs
mkdir -p $OUTDIR

DENOVO=../results/stacks/denovo
REFMAP=../results/stacks/refmap
FREEBA=../results/freebayes

vcftools --gzvcf $DENOVO/populations.snps.vcf --max-missing-count 116 --recode --out $OUTDIR/stacks_denovo --stdout | bgzip >$OUTDIR/stacks_denovo.vcf.gz
	vcftools --gzvcf $OUTDIR/stacks_denovo.vcf.gz --out $OUTDIR/stacks_denovo --missing-indv

vcftools --gzvcf $REFMAP/populations.snps.dict.vcf.gz --max-missing-count 116 --recode --out $OUTDIR/stacks_refmap --stdout | bgzip >$OUTDIR/stacks_refmap.vcf.gz
	vcftools --gzvcf $OUTDIR/stacks_refmap.vcf.gz --out $OUTDIR/stacks_refmap --missing-indv

vcfallelicprimitives $FREEBA/fb_parallel.vcf.gz | vcftools --vcf - --max-missing-count 116 --recode --out $OUTDIR/fb --stdout | bgzip >$OUTDIR/fb.vcf.gz
	vcftools --gzvcf $OUTDIR/fb.vcf.gz --out $OUTDIR/fb --missing-indv



# bcftools view $DENOVO/populations.snps.vcf | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/stacks_denovo.vcf.gz
# bcftools view $REFMAP/populations.snps.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/stacks_refmap.vcf.gz
# bcftools view $FREEBA/fb_parallel.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/fb.vcf.gz

# for file in $OUTDIR/*vcf.gz
# do tabix -f -p vcf $file
# done
