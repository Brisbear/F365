library(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(edgeR)
library(bedr)
library(tidyverse)
library(RColorBrewer)
library(plyr)
# library(GISTools)
library(tools)



# cytoband_df = cytoband$df
# chromosome = cytoband$chromosome
# xrange = c(cytoband$chr.len, 16569)
# mtband<-bed.df %>% filter(chr=='MT')

load("data_hk2/GRCh38.94.gtf.CDS.simple.RData")
treatment.list<-c('PMB','Colistin','F319','F287','F365')
load('degs.hk2.Rdata')
load('data_hk2/hg.genes.GRCh38.94.bed.RData')
load('data_hk2/crispr.genes.RData')
network<-read.table('data_hk2/signor_all_data_05_07_20.tsv',sep="\t",header=T,quote="",fill=T)
load('data_hk2/gene.alias.RData')
load('degs.cds.ann.RData')
genes<-degs.cds.ann[which(degs.cds.ann$Gene.Symbol %in% crisp.genes),c('Geneid','Gene.Symbol')]


out2<-cbind(out2,cpm(out2[,grep('^A\\d+$',colnames(out2))]))
colnames(out2)[(ncol(out2)-17):ncol(out2)]<-paste0(colnames(out2)[(ncol(out2)-17):ncol(out2)],'.cpm')
abundance<-apply(out2[rownames(degs[rowSums(degs)>=1,]),grep('^A\\d+.cpm$',colnames(out2))],1,mean)
abundance<-as.data.frame(abundance)


col_expr = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"),space='RGB')
a<-bed_new %>% filter(chr %in% c(1:22,'X','Y'))
a$chr<-paste0('chr',a$chr)
coordinate.degs<-a[which(a$ensembl_id %in% rownames(abundance)),c('chr','start','end','ensembl_id')]
abundance$ensembl<-rownames(abundance)
abundance<-left_join(coordinate.degs,abundance,by=c('ensembl_id'='ensembl'))
bed_new.ideogram<-a[,c('chr','start','end','ensembl_id')]
rownames(bed_new.ideogram)<-bed_new.ideogram$ensembl_id

# signor network
tg<-gene.alias %>% filter(gene.alias$Ensembl.Gene.ID %in% rownames(degs)[rowSums(degs)>=1])
tg<-gene.alias %>% filter(gene.alias$Ensembl.Gene.ID %in% rownames(degs)[rowSums(degs)>=1])
tg.source<- network %>% filter(network$ENTITYA %in% tg$Gene.Symbol)
tg.target<- network %>% filter(network$ENTITYB %in% tg$Gene.Symbol)
tg.net<-unique(rbind(tg.source,tg.target))
tg.net<-left_join(tg.net,gene.alias,by=c('ENTITYA'='Gene.Symbol'))
colnames(tg.net)[dim(tg.net)[2]]<-'ensembl_id.entityA'
tg.net<-left_join(tg.net,gene.alias,by=c('ENTITYB'='Gene.Symbol'))
colnames(tg.net)[dim(tg.net)[2]]<-'ensembl_id.entityB'
tg.net<-tg.net[!(is.na(tg.net$ensembl_id.entityA) | is.na(tg.net$ensembl_id.entityB)),]
tg.net.simple<-unique(tg.net[,c(9,10,29,30)])
tg.source.bed<- right_join(bed_new.ideogram, tg.net.simple,  by=c('ensembl_id'='ensembl_id.entityA'))
tg.target.bed<- right_join(bed_new.ideogram, tg.net.simple, by=c('ensembl_id'='ensembl_id.entityB'))
naindex<-is.na(tg.source.bed$chr) | is.na(tg.target.bed$chr)
tg.source.bed<-tg.source.bed[!naindex,]
tg.target.bed<-tg.target.bed[!naindex,]

fc<-merge(out2[rownames(degs[rowSums(degs)>=1,]),
      c(grep('(logFC)|(Geneid)',colnames(out2)))],bed_new.ideogram,by.x='Geneid',by.y='ensembl_id')
is.sig<-out2[which(out2$Geneid %in% fc$Geneid),grep('adj',colnames(out2))]<=0.05
expr<-fc[,2:6]
expr[!is.sig]<-NA
fc[,2:6]<-expr

fc<-fc[,c('chr','start','end','PMB.logFC', 'Colistin.logFC', 'F319.logFC', 'F287.logFC', 'F365.logFC')]


fc.index<-paste0(fc[,1],':',fc[,2],'-',fc[,3])
try(if (sum(is.valid.region(fc.index)) != dim(fc)[1]) stop('index is not valid...'))

fc.index.sort1 <- bedr.sort.region(fc.index,method='natural',engine = 'R',chr.to.num = c('X'=23,'Y'=24),check.chr = T)
fc<-fc[match(fc.index.sort1, fc.index),]


colnames(fc)[4:8]<-treatment.list

cytoband <- read.cytoband(species='hg38')
chromosomes<-cytoband$chromosome

# plot genomcis figure
pdf('circos.pdf',width=15,height=15)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species='hg38',plotType = NULL)

# text(0,0,'GRCh38',cex=1)
circos.genomicHeatmap(fc, col = col_expr, side = "outside",
                      border = NA,connection_height = convert_height(10, "mm"),
                      line_lwd=0.2,heatmap_height = 0.2)
circos.genomicIdeogram()

circos.genomicLabels(cbind(bed_new.ideogram[genes$Geneid,c(1:3)],genes$Gene.Symbol),
                     labels.column = 4, side = "inside",cex = 0.4)
circos.genomicDensity(abundance[,c(1:3,5)], col = c("#FF450080"), track.height = 0.1)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top", direction = "outside",labels.cex = 0.2*par("cex"))
})
mech<-unique(tg.source.bed$MECHANISM)
mm<-toTitleCase(as.character(mech))
mm[mm=='']<-'Unknown'

effect<-gsub('up.*','up',tg.source.bed$EFFECT)
effect<-gsub('down.*','down',effect)
effect<-mapvalues(as.factor(effect),from=c('up','down','unknown'),to=c('triangle','circle','ellipse'))

col_vector<-add.alpha(c(brewer.pal(9,'Set3'),brewer.pal(8,'Set1')),0.8)
tg.source.bed$MECHANISM<-as.factor(tg.source.bed$MECHANISM)


circos.genomicLink(tg.source.bed[,c(1:3)], tg.target.bed[,c(1:3)],lwd = 0.2*par("lwd"),
                   col=mapvalues(tg.source.bed$MECHANISM,from=mech,to=col_vector),
                   directional = 1, h.ratio=0.8,w=0.3,w2=0.1,
                   arr.length=0.1,arr.width = 0.05,border = NA)
                   # arr.type = effect)
circos.clear()

lgd_exp = Legend(at = seq(-5,5,1), col_fun = col_expr,
                   title_position = "topleft", title = "Log_2FC")
lgd_links = Legend(at = mm, type = "lines",
                    legend_gp = gpar(col = col_vector,lwd=3), title_position = "topleft",
                    title = "Interaction")
lgd_list_vertical = packLegend(lgd_exp, lgd_links)
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

dev.off()
