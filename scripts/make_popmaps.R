library(tidyverse)

# read in multiqc output
mqc <- read.table("../results/multiqc/multiqc_data/multiqc_fastqc.txt",sep="\t",header=TRUE)
mqc <- filter(mqc, str_detect(Sample,".1$")) %>% 
	mutate(., Sample=str_replace(Sample,regex(".1$"),"")) %>%
	select(Sample,Total.Sequences)
colnames(mqc) <-  gsub("\\.", "_",colnames(mqc))

# read in metadata table
meta <- read.table("../meta/Urban_ddRAD_FishIDs_Bioinformatics_2018.tsv",header=TRUE,sep="\t",quote="",comment.char="")

# add total # reads to metadata table, filter down to landscape genomics project and reads > 450k
meta <- left_join(x=meta,y=mqc,by=c("bioinformatics.id"="Sample")) %>% 
	filter(project=="LandscapeGenomics" & Total_Sequences > 450000)

# underscores only in column names
colnames(meta) <-  gsub("\\.", "_",colnames(meta))


# write out popmap files:

# first, a popmap file that subsamples for cstacks
	# we don't want to run all samples in cstacks, it takes too long
	# we'll take the 2 samples with the most reads from each population
	# this yields 56 samples

popmapcstacks <- group_by(meta, popmap_initial) %>% slice_max(Total_Sequences,n=2) %>% select(bioinformatics_id,popmap_initial)
write.table(popmapcstacks,"../meta/popmap_cstacks.txt", quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")

popmaptotal <- select(meta, bioinformatics_id,popmap_initial)
write.table(popmaptotal,"../meta/popmap_total.txt", quote=FALSE, col.names=FALSE,row.names=FALSE,sep="\t")