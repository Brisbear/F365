library(tidyverse)

load('GRCh38.94.gtf.RData')
bed_new<-df %>% separate(V9,c('GeneID','rest'),';')
bed_new$GeneID<-gsub('gene_id \\"(ENSG\\d+)\\"','\\1',bed_new$GeneID)

bed_new<-bed_new[,c(1:5,7,9)] %>% filter(V3=='gene')

colnames(bed_new)<-c('chr','gene_source','type','start','end','strand','ensembl_id')

save(bed_new,file='hg.genes.GRCh38.94.bed.RData')
