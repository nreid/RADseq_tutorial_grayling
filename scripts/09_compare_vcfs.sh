#!/bin/bash
#SBATCH --job-name=compare_vcfs
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


module load bcftools/1.9
module load htslib/1.9
module load vt/0.57721

INDIR=../results/filtered_vcfs

##############################
# check out filtering results:
##############################

vt peek $INDIR/fb_final.vcf.gz
vt peek $INDIR/denovo_final.vcf.gz
vt peek $INDIR/refmap_final.vcf.gz


##############################
# compare variant sets
##############################
# for two sets: VT partition
	# we can't directly compare the de novo set because it doesn't have genome coordinates
vt partition $INDIR/fb_final.vcf.gz $INDIR/refmap_final.vcf.gz


# extract intersections with bcftools isec
bcftools isec -O z -p $INDIR/isec $INDIR/fb_final.vcf.gz $INDIR/refmap_final.vcf.gz

# bcftools view -H isec_fb_gatk/0000.vcf | head -n 10

# look at a few incongruent markers from isec_fb_gatk/0000.vcf, which contains variants called by fb but not gatk

# chr20	31577045	.	T	C
	# called by fb, bcf, but not gatk. no clue why. 
	# apparently lots of support for C alternate allele. 
	# heterozygous parent has RO:35, AO:21


