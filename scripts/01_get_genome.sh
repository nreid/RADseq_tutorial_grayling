#!/bin/bash
#SBATCH --job-name=get_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# download genome assembly for Arctic grayling, Thymallus thymallus
# accession # GCA_004348285.1

GENOMEDIR=../genome
mkdir -p $GENOMEDIR

wget \
-P $GENOMEDIR \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/348/285/GCA_004348285.1_ASM434828v1/GCA_004348285.1_ASM434828v1_genomic.fna.gz


# index the genome using bwa

module load bwa/0.7.17
bwa index \
-p $GENOMEDIR/grayling \
$GENOMEDIR/GCA_004348285.1_ASM434828v1_genomic.fna.gz
