library(pheatmap)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(KEGGREST)

save_pheatmap_pdf <- function(x, filename, width=8, height=40) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

treatment<-c('PMB','Colistin','F319','F287','F365')
sample.list<-read.table('data_hk2/sample.final.tsv',sep="\t",header=T,quote="")
sample.list<-sample.list[order(factor(sample.list$Treatment,levels=treatment)),]
load('degs.hk2.Rdata')
load('data_hk2/gene.alias.RData')
load('data_hk2/ensembl.gene.name.GRch38.p13.RData')
load('data_hk2/GRCh38.94.gtf.CDS.simple.RData')
load('data_hk2/subloc.RData')
load('degs.cds.ann.RData')

# define gene set for  heatmap plot
# genes<-sort(c('ENSG00000077522','ENSG00000164400','ENSG00000143507','ENSG00000139318','ENSG00000065361','ENSG00000113578','ENSG00000066468','ENSG00000168621','ENSG00000113070','ENSG00000100385','ENSG00000171119','ENSG00000129946','ENSG00000160460'))
# genes<-degs.cds.ann$Geneid[grep('(^DNM)',degs.cds.ann$Gene.Symbol,perl=T)]
# genes<-degs.cds.ann$Geneid[grep('epsin',degs.cds.ann$gname,perl=T)]
# genes<-out2.cds$Geneid[which(out2.cds$F365.adj.P.Val<=0.05 & abs(out2.cds$F365.logFC)>=1)]

######################### Gene list from KEGG Pathway database #################
list.genes<-data.frame(unlist(keggLink("hsa",'path:hsa04144')))
list.genes<-gsub('hsa\\:','',list.genes[,1])
genes<-degs.cds.ann$Geneid[which(degs.cds.ann$entrez %in% list.genes)]

################### Gene list from target family ###############################
# load('target.family.RData')
# degs.cds.ann.family<-left_join(degs.cds.ann,target.family,by=c('Gene.Symbol'='HGNC.symbol'))
# a<-degs.cds.ann.family[which(rowSums(degs.cds.ann.family[,1:5])>=1 & !(is.na(degs.cds.ann.family$Type))),]
# a %>% group_by(Type) %>% summarize(n=n())
# genes<-a$Geneid[a$Type=='catalytic_receptor']

################### gene list from REACTOME ####################################
# table<-read.csv('any.cds.hk2.degs.reactome.csv')
# b<-unlist(strsplit(table[table$Pathway.identifier=='R-HSA-383280','Submitted.entities.found'], ";"))
# genes<-degs.cds.ann$Geneid[which(degs.cds.ann$entrez %in% b)]

################### Select genes based on toxicity order #######################
# genes<-out2.cds$Geneid[which(abs(out2.cds$PMB.logFC)>=abs(out2.cds$Colistin.logFC) &
#             abs(out2.cds$Colistin.logFC)>=abs(out2.cds$F319.logFC) &
#             abs(out2.cds$Colistin.logFC)>=abs(out2.cds$F287.logFC) &
#             abs(out2.cds$F319.logFC)>=abs(out2.cds$F365.logFC) &
#             abs(out2.cds$F287.logFC)>=abs(out2.cds$F365.logFC))]

# genes<-out2.cds$Geneid[which(abs(out2.cds$PMB.logFC)>=abs(out2.cds$Colistin.logFC) &
#             abs(out2.cds$Colistin.logFC)>=abs(out2.cds$F287.logFC) &
#             abs(out2.cds$F287.logFC)>=abs(out2.cds$F319.logFC) &
#             abs(out2.cds$F319.logFC)>=abs(out2.cds$F365.logFC))]

################## Select genes based on CRISPR results ########################
# load('crispr.genes.RData')
# print(crisp.genes)
# genes<-degs.cds.ann$Geneid[which(degs.cds.ann$Gene.Symbol %in% crisp.genes)]

item<-'Diff.Exp.Genes'
bkList = seq(-3, 3, by = 0.03)
expvalues<-out2.cds[which(out2.cds$Geneid %in% genes),grep('logFC',colnames(out2.cds))]

expvalues<-expvalues[order(rownames(expvalues)),]

y<-round(expvalues,2)
# degs.cds.loose<-data.frame(out2.cds[,grep('adj.P.Val',colnames(out2.cds))]<=0.05) & data.frame(abs(out2.cds[,grep('logFC',colnames(out2.cds))])>=log(1.5,2))
degs.cds.loose<-data.frame(out2.cds[,grep('adj.P.Val',colnames(out2.cds))]<=0.05)

degs.cds.loose<-1*degs.cds.loose
colnames(degs.cds.loose)<-treatment

# z<-degs.cds[which(out2.cds$Geneid %in% genes),c(1:5)]
z<-degs.cds.loose[which(out2.cds$Geneid %in% genes),c(1:5)]

z<-z[order(rownames(z)),]
y[z==0]<-''
y<-y[as.matrix(rowSums(z)>0),]

colnames(expvalues)<-gsub('.logFC','',colnames(expvalues))
g<-left_join(as.data.frame(rownames(expvalues)),gene.alias,by=c('rownames(expvalues)'='Ensembl.Gene.ID'))
g[is.na(g)]<-''
g<-g[as.matrix(rowSums(z)>0),]
gtmp<-left_join(g,gene.ensembl.ann,by=c('rownames(expvalues)'='ensembl'))
g1.rownames<-as.expression(lapply(as.matrix(unite(gtmp[,c(2,3,4)],'combined',sep=' / ')),function(x) bquote(italic(.(x)))))
expvalues<-expvalues[as.matrix(rowSums(z)>0),]
clustrow<-ifelse(nrow(expvalues)>1,T,F)

pp<-pheatmap(as.matrix(expvalues),
      color = colorRampPalette(c("blue","white","red"))(length(bkList)),
      cluster_cols=F,cluster_rows=clustrow,main=item,
      border_color='black',legend=T,labels_row = g1.rownames,
      show_rownames=T,cellwidth=8,cellheight=8,fontsize=4,
      breaks = bkList,display_numbers=y)

save_pheatmap_pdf(pp, paste0(item,'.pdf'))
