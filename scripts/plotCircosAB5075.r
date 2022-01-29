library(ComplexHeatmap)
library(circlize)
library(ape)
library(tidyverse)
library(RColorBrewer)

load('degs.5075.RData')
load('5075.rpkm.RData')
load('data_ab5075/genes.5075.RData')
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colnames(degs)<-paste0(colnames(degs),'.degs')
degs$Geneid<-rownames(degs)
degs$any.degs<-rowSums(degs[,1:16])


gff<-read.gff('data_ab5075/GCF_000963815.1_ASM96381v1_genomic.gff')
chr<-gff %>% filter(type=='region')
genes<-gff %>% filter(type=='gene')
genome<-data.frame(name=chr$seqid,start=chr$start,end=chr$end)
chr1<-chr %>% select(chr=seqid,start,end)
chr<- chr %>% select(seqid,end) %>% rename(length=end)

chr$add.len<-0

for (i in 2:nrow(genome)){
  genome$start[i]<-genome$end[i-1]+1
  genome$end[i]<-genome$end[i-1]+genome$end[i]
  chr$add.len[i]<-chr$add.len[i-1]+chr$length[i-1]
}
genes<-genes %>% separate(attributes,c('locus'),sep=";",extra = "drop") %>%
  mutate(locus = gsub('ID\\=gene-','',locus)) %>%
  left_join(chr,by=('seqid'))

expr<-out2 %>% select((contains('_Cont_') & contains('logFC')), Geneid) %>%
  filter (Geneid %in% genes$locus) %>%
  select(contains('L_') | contains('H_'), Geneid) %>%
  left_join(genes,by=c('Geneid'='locus')) %>%
  left_join(degs[,c(1:16,ncol(degs)-1,ncol(degs))],by=c('Geneid')) %>%
  filter(any.degs>=1)

  # mutate(start=start+add.len,end=end+add.len) %>%
expr1 <- expr %>% select(chr=seqid,start,end,contains('logFC'))

colnames(expr)<-gsub('.logFC','',colnames(expr))
gene.label<- expr %>% left_join(genes.5075,by=c('Geneid'='Accession.1')) %>% 
  select(chr=seqid,start,end,Gene.Name)
  # select(chr=seqid,start,end,Gene.Name)
# print(head(gene.label))
abundance<-data.frame(abundance=log10(apply(y_rpkm,1,mean)))
abundance$Geneid<-degs$Geneid
abundance<- abundance %>% filter(Geneid %in% degs[rowSums(degs[,1:16])>=1,'Geneid']) %>%
  left_join(expr,by=c('Geneid')) %>%
  select(chr=seqid,start,end,abundance)

col_expr = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"),space='RGB')
cog<-read.table('data_ab5075/GCF_000963815.1_ASM96381v1_protein.faa.eggnog.emapper.annotations.tsv',sep="\t")
cog<-cog %>% select(protein=V1,cog=V21)
feature<-read.table('data_ab5075/GCF_000963815.1_ASM96381v1_feature_table.txt',sep="\t",quote="")
feature <-feature %>% filter(V1=='CDS') %>% select(gene=V17,protein=V11) %>%
  left_join(cog,by=c('protein')) %>% 
  left_join(genes,by=c('gene'='locus')) %>%
  select(chr=seqid,start,end,cog) %>%
  mutate(cog=replace(cog,is.na(cog),'')) %>%
  mutate(cog=str_sub(cog,1,1))

color.cog<-data.frame(cog=sort(unique(feature$cog)), 
                      color=getPalette(length(unique(feature$cog))))
color.cog<-feature %>% left_join(color.cog,by=c('cog'))



# circlos plot
pdf('circos.5075.pdf',width=15,height=15)
circos.par(start.degree = 90,track.margin=c(0,0))
# circos.genomicInitialize(genome)
circos.genomicInitialize(chr1)
circos.genomicHeatmap(expr1, col = col_expr, side = "outside",
                      border = NA,connection_height = convert_height(10, "mm"),
                      line_lwd=0.05,heatmap_height = 0.2)
# circos.genomicTrackPlotRegion(color.cog[,c(1:3,5)],ylim=c(0,1),
#                     panel.fun = function(region,value,...){
#                       circos.genomicRect(region,value,ytop.column = 1,ybottom = 0,
#                                          col==value )})
circos.genomicLabels(gene.label,labels.column = 4, side = "inside",cex = 0.05,line_lwd=0.05)
circos.genomicDensity(abundance, col = c("#E3242B80"), track.height = 0.2)
circos.clear()
dev.off()




