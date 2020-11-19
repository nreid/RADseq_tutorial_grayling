library(tidyverse)

# this script aggregates 
	# 1) metadata
	# 2) mapping data
	# 3) barcode/pool data
	

# read in metadata table
	meta <- read.table("../meta/Urban_ddRAD_FishIDs_Bioinformatics_2018.tsv",header=TRUE,sep="\t",quote="",comment.char="")

# read in sam stats table

	sn <- read.table("../results/align_stats/SN.txt",sep="\t") %>% t()
	sn <- cbind(Sample=rownames(sn),data.frame(sn))
	colnames(sn) <- colnames(sn) %>% str_replace(regex("\\.*$"),"")

	sn <- sn[,c(1,2,8,9,14,20,24,39)]

# add mapping data to metadata table
	meta <- left_join(x=meta,y=sn,by=c("bioinformatics.id"="Sample")) 

# underscores only in column names
	colnames(meta) <-  gsub("\\.", "_",colnames(meta))


# read in barcode and pool information
	f <- list.files("../meta",pattern="barcode*",full.names=TRUE)

	bc <- data.frame()
	for(i in f){
		seqpool <- str_extract(i,regex("[0-9][0-9][0-9]"))
		print(pool)
		bc <- rbind(bc,cbind(pool=seqpool,read.table(i)))
	}

	colnames(bc) <- c("pool","barcode","Sample")

# merge barcode and pool info with metadata and mapping data

	meta <- left_join(x=meta,y=bc,by=c("bioinformatics_id"="Sample")) 

# have a look

# a bunch of individuals in pools 1 and 2 have major problems. these will get dropped when filtering the VCF
plot(meta$reads_mapped/meta$raw_total_sequences,col=factor(meta$pool)

