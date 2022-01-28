library(pheatmap)
# library(stringr)
# library(dplyr)
library(tidyverse)
library(ggplot2)
library(KEGGREST)
library(edgeR)

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

out2<-cbind(out2,cpm(out2[,grep('^A\\d+$',colnames(out2))]))
colnames(out2)[(ncol(out2)-17):ncol(out2)]<-paste0(colnames(out2)[(ncol(out2)-17):ncol(out2)],'.cpm')
out2.cds<-cbind(out2.cds,out2[rownames(out2.cds),(ncol(out2)-17):ncol(out2)])
load('data_hk2/gene.alias.RData')
load('data_hk2/ensembl.gene.name.GRch38.p13.RData')
load('data_hk2/GRCh38.94.gtf.CDS.simple.RData')
load('data_hk2/subloc.RData')
load('degs.cds.ann.RData')

# genes<-degs.cds.ann$Geneid[grep('(^ATG)',degs.cds.ann$Gene.Symbol,perl=T)]
genes<-degs.cds.ann$Geneid[grep('potassium',degs.cds.ann$gname,perl=T)]

# list.genes<-data.frame(unlist(keggLink("hsa",'path:hsa04144')))
# list.genes<-gsub('hsa\\:','',list.genes[,1])
# genes<-degs.cds.ann$Geneid[which(degs.cds.ann$entrez %in% list.genes)]

# load('target.family.RData')
# degs.cds.ann.family<-left_join(degs.cds.ann,target.family,by=c('Gene.Symbol'='HGNC.symbol'))
# a<-degs.cds.ann.family[which(rowSums(degs.cds.ann.family[,1:5])>=1 & !(is.na(degs.cds.ann.family$Type))),]
# a %>% group_by(Type) %>% summarize(n=n())
# genes<-a$Geneid[a$Type=='catalytic_receptor']

# table<-read.csv('any.cds.hk2.degs.reactome.csv')
# b<-unlist(strsplit(table[table$Pathway.identifier=='R-HSA-383280','Submitted.entities.found'], ";"))
# genes<-degs.cds.ann$Geneid[which(degs.cds.ann$entrez %in% b)]

# genes<-out2.cds$Geneid[which(abs(out2.cds$PMB.logFC)>=abs(out2.cds$Colistin.logFC) &
#             abs(out2.cds$Colistin.logFC)>=abs(out2.cds$F319.logFC) &
#             abs(out2.cds$Colistin.logFC)>=abs(out2.cds$F287.logFC) &
#             abs(out2.cds$F319.logFC)>=abs(out2.cds$F365.logFC) &
#             abs(out2.cds$F287.logFC)>=abs(out2.cds$F365.logFC))]

# genes<-out2.cds$Geneid[which(abs(out2.cds$PMB.logFC)>=abs(out2.cds$Colistin.logFC) &
#             abs(out2.cds$Colistin.logFC)>=abs(out2.cds$F287.logFC) &
#             abs(out2.cds$F287.logFC)>=abs(out2.cds$F319.logFC) &
#             abs(out2.cds$F319.logFC)>=abs(out2.cds$F365.logFC))]
# genes<-out2.cds$Geneid[which(out2.cds$F365.adj.P.Val<=0.05 & abs(out2.cds$F365.logFC)>=1)]

# load('crispr.genes.RData')
# print(crisp.genes)
# genes<-degs.cds.ann$Geneid[which(degs.cds.ann$Gene.Symbol %in% crisp.genes)]

item<-'CPM'
item1<-'Read.counts'
expvalues<-out2.cds[which(out2.cds$Geneid %in% genes),grep('^A\\d+.cpm$',colnames(out2.cds))]
expvalues1<-out2.cds[which(out2.cds$Geneid %in% genes),grep('^A\\d+$',colnames(out2.cds))]

coln<-colnames(expvalues1)
expvalues<-expvalues[,paste0(coln[match(sample.list$Sample.ID,coln)],'.cpm')]
expvalues1<-expvalues1[,coln[match(sample.list$Sample.ID,coln)]]

expvalues<-expvalues[order(rownames(expvalues)),]
expvalues1<-expvalues1[order(rownames(expvalues1)),]

y<-round(expvalues,2)
y1<-expvalues1
# degs.cds.loose<-data.frame(out2.cds[,grep('adj.P.Val',colnames(out2.cds))]<=0.05) & data.frame(abs(out2.cds[,grep('logFC',colnames(out2.cds))])>=log(1.5,2))
degs.cds.loose<-data.frame(out2.cds[,grep('adj.P.Val',colnames(out2.cds))]<=0.05)
degs.cds.loose<-1*degs.cds.loose
colnames(degs.cds.loose)<-treatment

# z<-degs.cds[which(out2.cds$Geneid %in% genes),c(1:5)]
z<-degs.cds.loose[which(out2.cds$Geneid %in% genes),c(1:5)]

z<-z[order(rownames(z)),]
print(dim(z))
print(head(z))

y<-y[as.matrix(rowSums(z)>0),]
y1<-y1[as.matrix(rowSums(z)>0),]

g<-left_join(as.data.frame(rownames(expvalues)),gene.alias,by=c('rownames(expvalues)'='Ensembl.Gene.ID'))
g[is.na(g)]<-''
g<-g[as.matrix(rowSums(z)>0),]
gtmp<-left_join(g,gene.ensembl.ann,by=c('rownames(expvalues)'='ensembl'))
g1.rownames<-as.expression(lapply(as.matrix(unite(gtmp[,c(2,3,4)],'combined',sep=' / ')),function(x) bquote(italic(.(x)))))
expvalues<-expvalues[as.matrix(rowSums(z)>0),]
expvalues1<-expvalues1[as.matrix(rowSums(z)>0),]
bkList = seq(0, max(expvalues), by = max(expvalues)/200)
bkList1 = seq(0, max(expvalues1), by = max(expvalues)/200)
clustrow<-ifelse(nrow(expvalues)>1,T,F)

annotation_col<-data.frame(Condition=sample.list$Treatment)
rownames(annotation_col)<-paste0(sample.list$Sample.ID,'.cpm')
p1<-pheatmap(as.matrix(expvalues),annotation_col = annotation_col,
      color = colorRampPalette(c("white","purple"))(length(bkList)),cluster_cols=F,cluster_rows=F,
      border_color='black',legend=T,
      plotScale='log2',
      labels_row = g1.rownames,main=item,
      show_rownames=T,cellwidth=8,cellheight=8,fontsize=6,breaks = bkList,
      display_numbers=y)
save_pheatmap_pdf(p1, paste0(item,'.pdf'))
dev.off()
p2<-pheatmap(as.matrix(expvalues1),annotation_col = annotation_col,
      color = colorRampPalette(c("white","blue"))(length(bkList1)),cluster_cols=F,cluster_rows=clustrow,
      border_color='black',legend=T,
      plotScale='log2',
      labels_row = g1.rownames,main=item1,
      show_rownames=T,cellwidth=8,cellheight=8,fontsize=4,breaks = bkList1,
      display_numbers=y1)

save_pheatmap_pdf(p2, paste0(item1,'.pdf'))
