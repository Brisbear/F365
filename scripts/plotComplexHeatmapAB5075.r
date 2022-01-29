# library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(readxl)



load('degs.5075.RData')
load('5075.rpkm.RData')
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
y_rpkm<-as.data.frame(rowMeans(y_rpkm))
colnames(y_rpkm)<-'rpkm'

genes<-read.table('data_ab5075/All_genes_of_A._baumannii_AB5075-UW.txt',header=T,sep="\t",fill=T,quote="")
bkList = seq(-5, 5, by = 0.05)
rownames(out2)<-out2$Geneid

biocyc.logFC<-read.table('combined_logFC_ab.txt',header=T,sep="\t") 

df0<-read_excel('data_ab5075/selected_genes_for_heatmap.xlsx',sheet='Sheet1') %>%
  filter(grepl('(Fatty Acid)|(lipid A)|(Aromatic)|(Respiration)',pathway)) %>% 
  left_join(genes,by=c('gene'='Accession.1'))

df1<-df0 %>% 
  left_join((out2 %>% dplyr::select(colnames(biocyc.logFC)[2:17],Geneid)),by=c('gene'='Geneid')) %>%
  dplyr::select(contains('L_') | contains('H_'), gene,Gene.Name, Product, pathway) 
df1$pathway<-as.factor(df1$pathway)


df2<-as.matrix(round(df1[,1:16],2))
a<-df1 %>%
  left_join(degs %>% mutate(Geneid = rownames(degs)), by=c('gene'='Geneid')) %>%
  dplyr::select((contains('L_') | contains('H_')) & ! contains('logFC')  & contains('_Cont_'))

df2[a==0]<-''

make_italics <- function(x) {
  if (grepl('ABUW',x) ==T){
    bquote(.(x))
  }else{
    bquote(italic(.(x)))
  }
}

b<-df1 %>%
  dplyr::select(gene,Gene.Name)

b1<-as.expression(lapply(b$Gene.Name,function(x) make_italics(x)))

ha = rowAnnotation(FPKM = anno_barplot(log10(y_rpkm[b$gene,'rpkm'])), annotation_name_rot = 90)
colourCount <- length(unique(df1$pathway))
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
pdf(paste0("selected.expression.pdf"), width = 10, height = 25)
ht1<-Heatmap(as.matrix(df1[,1:16]), name = "log2FC", col = col_fun,
             rect_gp = gpar(col = "black", lwd = 0.5),
             column_title = "Treatment",column_title_gp = gpar(fontsize = 10, fontface = "bold"),
             row_title = "Gene",row_title_gp = gpar(fontsize = 11, fontface = "bold"),
             cluster_columns = F,
             left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill=getPalette(colourCount)),
                                                              labels = LETTERS[1:length(unique(df1$pathway))],
                                                               labels_gp = gpar(col = "black", fontsize = 10))),
             # show_row_dend =F,
             row_names_gp = gpar(fontsize = 4),row_labels = b1,
             column_names_gp=gpar(fontsize=6),column_labels = gsub('^(\\w+).(\\w)_Cont_(\\w)h.logFC','\\1-\\2-\\3h',colnames(df1)[1:16]),
             right_annotation = ha,
             width = 16*unit(2, "mm"),
             height= nrow(df1)*unit(2, "mm")
             # cell_fun = function(j, i, x, y, w, h, col) {
             #   grid.text(sprintf('%s',df2[i, j]), x, y,  gp = gpar(fontsize = 1.8)
             #   )}
             )
draw(ht1,row_split = df1$pathway,ht_gap=unit(2,'mm'))
dev.off()
