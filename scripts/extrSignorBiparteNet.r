library(tidyverse)

load('degs.hk2.Rdata')
load('data_hk2/hg.genes.GRCh38.94.bed.RData')
a<-bed_new %>% filter(chr %in% c(1:22,'X','Y'))
a$chr<-paste0('chr',a$chr)
bed_new.ideogram<-a[,c('chr','start','end','ensembl_id')]

network<-read.table('data_hk2/signor_all_data_05_07_20.tsv',sep="\t",header=T,quote="",fill=T)
load('data_hk2/gene.alias.RData')

tg<-gene.alias %>% filter(gene.alias$Ensembl.Gene.ID %in% rownames(degs)[rowSums(degs)>=1])

# tg.source<- tg %>% filter(tg$Gene.Symbol %in% network$ENTITYA)
tg.source<- network %>% filter(network$ENTITYA %in% tg$Gene.Symbol)
# tg.target<- tg %>% filter(tg$Gene.Symbol %in% network$ENTITYB)
tg.target<- network %>% filter(network$ENTITYB %in% tg$Gene.Symbol)

tg.net<-unique(rbind(tg.source,tg.target))
tg.net<-left_join(tg.net,gene.alias,by=c('ENTITYA'='Gene.Symbol'))
colnames(tg.net)[dim(tg.net)[2]]<-'ensembl_id.entityA'
tg.net<-left_join(tg.net,gene.alias,by=c('ENTITYB'='Gene.Symbol'))
colnames(tg.net)[dim(tg.net)[2]]<-'ensembl_id.entityB'
tg.net<-tg.net[!(is.na(tg.net$ensembl_id.entityA) | is.na(tg.net$ensembl_id.entityB)),]

tg.source.bed<- right_join(bed_new.ideogram, tg.net, by=c('ensembl_id'='ensembl_id.entityA'))
tg.target.bed<- right_join(bed_new.ideogram, tg.net, by=c('ensembl_id'='ensembl_id.entityB'))
