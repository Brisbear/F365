library(ComplexHeatmap)
library(stringr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(KEGGREST)
library(circlize)
library(dendsort)
library(edgeR)
# library(randomcoloR)
library(RColorBrewer)

treatment<-c('PMB','Colistin','F319','F287','F365')
sample.list<-read.table('data_hk2/sample.final.tsv',sep="\t",header=T,quote="")
sample.list<-sample.list[order(factor(sample.list$Treatment,levels=treatment)),]
load('degs.hk2.Rdata')
load('data_hk2/gene.alias.RData')
load('data_hk2/ensembl.gene.name.GRch38.p13.RData')
load('data_hk2/GRCh38.94.gtf.CDS.simple.RData')
load('data_hk2/hk2.rpkm.RData')
load('data_hk2/transporterClass.RData')


######################### Gene list from KEGG Pathway database ################
# list.genes<-data.frame(unlist(keggLink("hsa",'path:hsa04210')))
# list.genes<-gsub('hsa\\:','',list.genes[,1])
# genes<-degs.cds.ann$Geneid[which(degs.cds.ann$entrez %in% list.genes)]
# df<-read.table('ubiquitination/ubiquitination_genes.ensgonly.txt',header=F)
# genes<-unique(df$V1)
# genes<-degs.cds.ann$Geneid[grep('^ZFP',degs.cds.ann$Gene.Symbol,perl=T)]
# keltch<-read.table('ubiquitination/keltch_containing_prot.txt',header=F)
# genes<-degs.cds.ann[which(degs.cds.ann$Gene.Symbol %in% keltch$V1), 'Geneid']
# ddb2<-read.table('ubiquitination/exported-data.txt',sep="\t",header=T)
# genes<-degs.cds.ann[which(degs.cds.ann$Gene.Symbol %in% ddb2[ddb2$SCORE>0.65,'SUBGENE']), 'Geneid']
# genes<-degs.cds.ann$Geneid[grep('SLC',degs.cds.ann$gname,perl=T)]
# load('apoptosis/apoptosis.genes.RData')
# apop<-apoptosis_genes[which(apoptosis_genes$species=='HUMAN'),'external_id']
# genes<-degs.cds.ann[which(degs.cds.ann$Gene.Symbol %in% apop), 'Geneid']

geneinfo<-read_excel('data_hk2/geneinfor.xlsx',sheet=1) %>%
  mutate(class=ifelse(grepl('Channel|channel',class),'Ion channel',class))
geneinfo<-geneinfo[!duplicated(geneinfo$ensembl),]
genes<-geneinfo$ensembl

degs.cds.loose<-data.frame(out2.cds[,grep('adj.P.Val',colnames(out2.cds))]<=0.05) &
  data.frame(abs(out2.cds[,grep('logFC',colnames(out2.cds))])>=log(1.5,2))
# degs.cds.loose<-data.frame(out2.cds[,grep('adj.P.Val',colnames(out2.cds))]<=0.05)
degs.cds.loose<-1*degs.cds.loose
colnames(degs.cds.loose)<-treatment

# sig.trpts<-trpt %>% filter(ensembl %in% rownames(degs.cds.loose)[rowSums(degs.cds.loose)>=1])
# sig.trpts <-filter(sig.trpts,grepl('channel|SLC|ATPase|transporter',Family.name))

expvalues<-out2.cds[which(out2.cds$Geneid %in% genes),grep('logFC',colnames(out2.cds))]
expvalues<-expvalues[order(rownames(expvalues)),]
abundance<-y_rpkm[which(rownames(y_rpkm) %in% genes),]
abundance<-abundance[order(rownames(abundance)),match(sample.list$Sample.ID,colnames(abundance))]

degs.cds.loose[which(rownames(degs.cds.loose) %in% rownames(abundance)[apply(abundance,1,max)<1]),]<-0
# z<-degs.cds[which(out2.cds$Geneid %in% genes),c(1:5)]
z<-degs.cds.loose[which(out2.cds$Geneid %in% genes),c(1:5)]
z<-z[order(rownames(z)),]
abundance<-abundance[as.matrix(rowSums(z)>0),]

colnames(expvalues)<-gsub('.logFC','',colnames(expvalues))
num.text<-round(expvalues,2)
num.text[z==0]<-''
num.text<-num.text[as.matrix(rowSums(z)>0),]
g<-left_join(as.data.frame(rownames(expvalues)),gene.alias,by=c('rownames(expvalues)'='Ensembl.Gene.ID'))
g[is.na(g)]<-''
g<-g[as.matrix(rowSums(z)>0),]
gtmp<-left_join(g,gene.ensembl.ann,by=c('rownames(expvalues)'='ensembl'))
# g1.rownames<-as.expression(lapply(as.matrix(unite(gtmp[,c(2,3,4)],'combined',sep=' / ')),function(x) bquote(italic(.(x)))))
g1.rownames<-as.expression(lapply(as.matrix(gtmp[,2]),function(x) bquote(italic(.(x)))))

expvalues<-expvalues[rowSums(z)>0,]
clustrow<-ifelse(nrow(expvalues)>1,T,F)

load('data_hk2/subloc_df.RData')
sl<-as.matrix(d[rownames(expvalues),])
rownames(sl)<-rownames(expvalues)
sl[is.na(sl)] <- 0
sl<-sl[,colSums(sl)!=0]

# network<-read.table('data_hk2/signor_all_data_06_09_20.tsv',sep="\t",header=T,quote="",fill=T)
# save(network,file='signor_net.0906.RData')
load('data_hk2/signor_net.0906.RData')
degr<-data.frame(degree=unlist(lapply(gtmp$Gene.Symbol, function(x) length(network[which(network$ENTITYA %in% x | network$ENTITYB %in% x),1]))))
rownames(degr)<-gtmp[,'rownames(expvalues)']

col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_fun_anno = colorRamp2(c(0, 200), c("white", "darkred"))
dend = dendsort(hclust(dist(as.matrix(expvalues))))
avg.cpm<-rowSums(abundance)
pdf("heatmapSelectGenes_HK2.pdf", width = 10, height = 20)
#
# uucd<-read.table('ubiquitination/classification_uucd.txt',sep="\t",header=T,quote="",fill=T)
# ann.uucd<-lapply(uucd$Annotation,function(x) str_split(x, '; '))
# for (i in 1:length(ann.uucd)){
#   ann.uucd[i]<-unlist(ann.uucd[[i]])[1]
# }
# uucd$Annotation<-unlist(ann.uucd)
# gtmp<-left_join(gtmp,uucd,by=c('Gene.Symbol'='Annotation'))
gtmp[is.na(gtmp)]<-' '

# palette <- data.frame(family=unique(gtmp$family),color=distinctColorPalette(length(unique(gtmp$family))))
# uucd.col<-left_join(gtmp,palette,by=c('family'='family'))

# with UUCD
# ha = HeatmapAnnotation(Degree=degr[,1],UUCD=gtmp$family,Abundance = anno_barplot(apply(abundance,1,mean)),
#     annotation_name_rot = c(90, 90,90),col = list(Degree = col_fun_anno,UUCD=setNames(uucd.col$color,gtmp$family)),which='row')
# without UUCD
a<-left_join(data.frame(name=rownames(expvalues)),geneinfo,by=c('name'='ensembl'))
a$class<-as.factor(a$class)
logscale.breaks<-c(1,10,100,1000,10000)
ha = HeatmapAnnotation(Degree=degr[,1],RPKM = anno_barplot(log10(apply(abundance,1,mean)),
                                                                 axis_param = list(at = log10(logscale.breaks), labels =logscale.breaks,gp = gpar(fontsize=3))),
    annotation_name_rot = c(90, 90,90),col = list(Degree = col_fun_anno),which='row')

colourCount <- length(unique(a$class))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ht1<-Heatmap(as.matrix(expvalues), name = "log2FC", col = col_fun,
        rect_gp = gpar(col = "black", lwd = 0.5),

        column_title = "Treatment",column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        row_title = "Gene",row_title_gp = gpar(fontsize = 20, fontface = "bold"),
        cluster_columns = F,
        # cluster_rows = dend,
        row_names_gp = gpar(fontsize = 6),row_labels = g1.rownames,
        column_names_gp=gpar(fontsize=6),
        right_annotation = ha,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill=getPalette(colourCount)),
                                                         labels = LETTERS[1:length(unique(a$class))],
                                                         labels_gp = gpar(col = "black", fontsize = 10))),
       #  cell_fun = function(j, i, x, y, w, h, col) {
       #  grid.text(sprintf('%s',num.text[i, j]), x, y,  gp = gpar(fontsize = 4)
       # )},
        width = 5*unit(3, "mm"),
       height = nrow(expvalues)*unit(3,'mm'))
a.form<-data.frame(id=LETTERS[1:length(unique(a$class))],class=unique(a$class))
save(a.form,file='gene.class.heatmap.RData')
ht2 = Heatmap(sl, name = "subcellular localisation",col=c('white','darkgreen'),
        column_title = "Location",column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        row_title = "Gene",row_title_gp = gpar(fontsize = 20, fontface = "bold"),
        row_order = dend$order,cluster_columns = F,row_labels = g1.rownames,
        column_names_gp=gpar(fontsize=6),show_heatmap_legend=F,
        row_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "black", lwd = 1),
        width = unit(2, "cm"))
# ht3 = Heatmap(tc, name = "Levels in kidney tubules", col=structure(1:5,names=c("Unknown","Not detected",'Low','Medium','High')),
#         column_title = "Levels in kidney tubules",column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#         row_title = "Gene",row_title_gp = gpar(fontsize = 20, fontface = "bold"),
#         row_order = dend$order,cluster_columns = F,row_labels = g1.rownames,
#         column_names_gp=gpar(fontsize=6),show_heatmap_legend=F,
#         row_names_gp = gpar(fontsize = 6),
#         rect_gp = gpar(col = "black", lwd = 1),
#         width = unit(2, "cm"))

htlist<-ht1
draw(htlist,row_split = a$class,ht_gap=unit(5,'mm'))
dev.off()
