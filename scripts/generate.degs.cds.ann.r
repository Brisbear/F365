library(tidyverse)


load('data_hk2/ensembl.gene.name.GRch38.p13.RData')
load('data_hk2/gene.alias.RData')
genes.cds.ann<-gene.alias %>% full_join(gene.ensembl.ann,by=c('Ensembl.Gene.ID'='ensembl'))

load('degs.hk2.RData')
degs.cds.ann<-degs %>% rownames_to_column('Geneid') %>%
	left_join(genes.cds.ann,by=c('Geneid'='Ensembl.Gene.ID'))

save(degs.cds.ann,file='data_hk2/degs.cds.ann.RData')
