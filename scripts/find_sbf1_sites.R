library(tidyverse)
library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()


# This R script finds sbf1 cut sites, flanking sites and writes them to bed files

# sbf1 <- "CCTGCAGG"
# msei <- "TTAA"

# read in reference genome 
seqs <- readDNAStringSet("../genome/GCA_004348285.1_ASM434828v1_genomic.fna")

# find sbf1 sites
mp <- vmatchPattern("CCTGCAGG",seqs)

# make it a GRanges object
mp2 <- as(mp, "GRanges")

# make it a data frame
mp3 <- as.data.frame(mp2)
mp3[,1] <- str_replace(mp3[,1],regex(" .*"),"")

# 0-index start
mp3[,2] <- mp3[,2] - 1

# left flank
left <- mp3
	left[,2] <- left[,2] - 9
	left[,3] <- left[,3] - 9
left[left[,2] < 0,2] <- 0

# right flank
right <- mp3
	right[,2] <- right[,2] + 9
	right[,3] <- right[,3] + 9

# find sbf1 1-off sites
mp_1 <- vmatchPattern("CCTGCAGG",seqs,max.mismatch=1)

# make it a GRanges object
mp2_1 <- as(mp_1, "GRanges")

# make it a data frame
mp3_1 <- as.data.frame(mp2_1)
mp3_1[,1] <- str_replace(mp3_1[,1],regex(" .*"),"")

# 0-index start
mp3_1[,2] <- mp3_1[,2] - 1


# write tables
write.table(mp3[,1:3],"../meta/sbf1.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(mp3_1[,1:3],"../meta/sbf1_off.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)


# find the distances from sbf1 to mseI cut sites

# mseI: TTAA

# find msei sites and get granges object
msei <- vmatchPattern("TTAA",seqs)
msei2 <- as(msei, "GRanges")

# get closest upstream and downstream msei cut sites to sbf1 cut sites
rightind <- precede(mp2,msei2,"first")
	NOrightNA <- !is.na(rightind)
leftind <- follow(mp2,msei2,"last")
	NOleftNA <- !is.na(leftind)

out <- matrix(nrow=length(mp2),ncol=2)

out[NOrightNA,1] <- distance(mp2[NOrightNA,],msei2[rightind[NOrightNA],])
out[NOleftNA,2] <- distance(mp2[NOleftNA,],msei2[leftind[NOleftNA],])

# `out` now contains all SbfI-MseI RAD fragment lengths in the reference genome. 
	# column 1 is the left side of SbfI, column 2 is the right. 
	
