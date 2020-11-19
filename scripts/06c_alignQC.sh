#!/bin/bash 
#SBATCH --job-name=samstats
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --qos=general
#SBATCH --partition=general

# calculate alignment statistic for each file

module load samtools/1.10

# input, output directories

INDIR=../results/aligned
OUTDIR=../results/align_stats
mkdir -p $OUTDIR


# samtools bam statistics
for file in $(find $INDIR -name "*bam"); 
do 
SAM=$(basename $file .bam)
samtools stats $file >$OUTDIR/${SAM}.stats
echo $file;
done


# put the basic stats all in one file. 

FILES=($(find $OUTDIR -name "*.stats" | sort))

grep "^SN" ${FILES[0]} | cut -f 2 > $OUTDIR/SN.txt
for file in ${FILES[@]}
do paste $OUTDIR/SN.txt <(grep ^SN $file | cut -f 3) > $OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt
	echo $file
done

# add a header with sample names
cat \
<(echo ${FILES[@]} | sed 's,../results/align_stats/,,g' | sed 's/.stats//g' | sed 's/ /\t/g') \
$OUTDIR/SN.txt \
>$OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt
