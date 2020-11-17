
#!/bin/bash
#SBATCH --job-name=raw_data_symlinks
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=128M
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# this script symlinks the raw data pools to the working directory

RAWDATADIR=../data/pools
mkdir -p $RAWDATADIR

POOLPATH="/labs/Wegrzyn/Urban_RAD_ArcticGrayling/raw_fastq/raw_fastq_files"

for f in ${fpath}/Golden-Pool*fastq.gz; do
        echo $f
        echo `basename ${f}`
        ln -s ${f} `basename ${f}`
done

