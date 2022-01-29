library(limma)
library(tidyverse)
library(edgeR)
library(jsonlite)
library(mixOmics)
library(RColorBrewer)
library(pheatmap)
library(UpSetR)
library(ComplexHeatmap)


cols <- brewer.pal(8, "Paired")
cols<-append(cols, c("#808080"), 2)

counts_file <- 'data_ab5075/merged.ab.counts.txt'
s<-read.table('data_ab5075/sample.tsv',sep = "\t",header=T)
s$Sample.ID<-gsub('_L.*','',s$Sample.ID)
count_cols <- s$Sample.ID
s$Group<-paste(s$Antibiotic,ifelse(s$Concentration==2,'L','H'),paste0(s$Time,'h'),sep="_")
s$Group[grep('Control',s$Group)]<-gsub('_H_','_',s$Group[grep('Control',s$Group)])
# conditions
g<-c('Control_0h','Control_1h','Colistin_L_1h',
  'Colistin_H_1h','F287_L_1h','F287_H_1h','F319_L_1h','F319_H_1h',
  'F365_L_1h','F365_H_1h','Control_4h','Colistin_L_4h','Colistin_H_4h',
  'F287_L_4h','F287_H_4h','F319_L_4h','F319_H_4h','F365_L_4h','F365_H_4h')
# comparisons
treatment.list<-c('Col.L_Cont_1h','Col.H_Cont_1h','F287.L_Cont_1h','F287.H_Cont_1h',
'F319.L_Cont_1h','F319.H_Cont_1h','F365.L_Cont_1h','F365.H_Cont_1h','Col.L_Cont_4h',
'Col.H_Cont_4h','F287.L_Cont_4h','F287.H_Cont_4h','F319.L_Cont_4h','F319.H_Cont_4h','F365.L_Cont_4h','F365.H_Cont_4h',
'Cont_1h_0h','Cont_4h_0h','Cont_4h_1h',
'Col.L_4h_1h','Col.H_4h_1h','F287.L_4h_1h','F287.H_4h_1h',
'F319.L_4h_1h','F319.H_4h_1h','F365.L_4h_1h','F365.H_4h_1h')
design<-matrix(c(as.numeric(grepl('Control_0h',s$Group)),
                  as.numeric(grepl('Control_1h',s$Group)),
                  as.numeric(grepl('Colistin_L_1h',s$Group)),
                  as.numeric(grepl('Colistin_H_1h',s$Group)),
                  as.numeric(grepl('F287_L_1h',s$Group)),
                  as.numeric(grepl('F287_H_1h',s$Group)),
                  as.numeric(grepl('F319_L_1h',s$Group)),
                  as.numeric(grepl('F319_H_1h',s$Group)),
                  as.numeric(grepl('F365_L_1h',s$Group)),
                  as.numeric(grepl('F365_H_1h',s$Group)),
                  as.numeric(grepl('Control_4h',s$Group)),
                  as.numeric(grepl('Colistin_L_4h',s$Group)),
                  as.numeric(grepl('Colistin_H_4h',s$Group)),
                  as.numeric(grepl('F287_L_4h',s$Group)),
                  as.numeric(grepl('F287_H_4h',s$Group)),
                  as.numeric(grepl('F319_L_4h',s$Group)),
                  as.numeric(grepl('F319_H_4h',s$Group)),
                  as.numeric(grepl('F365_L_4h',s$Group)),
                  as.numeric(grepl('F365_H_4h',s$Group))),ncol=19,
dimnames=list(s$Sample.ID,g))

cont.matrix <- matrix(c(c(0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0),
                        c(0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,1),
                        c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                        c(-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                        c(0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                        c(0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
                        c(0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                        c(0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
                        c(0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0),
                        c(0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0),
                        c(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0),
                        c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1)
                      ),ncol=27,dimnames=list(g,treatment.list))

export_cols <- c(c('Geneid','Chr','Start','End','Strand','Length'),s$Sample.ID[order(s$Group)])

x<-read.delim(counts_file,
              sep="	",
              check.names=FALSE,
              colClasses='character',
              na.strings=c(),
              skip=0
              )

# Now re-read the first header line.  Workaround R problem that always has strip.white=T for read.table
colnames(x) <- scan(counts_file,
                    what="",
                    sep="	",
                    nlines=1,
                    strip.white=F,
                    quote = "\"",
                    skip="0"
                    )

x[,count_cols] <- apply(x[,count_cols], 2, function(v) as.numeric(v))     # Force numeric count columns
counts <- x[, count_cols]


keepMin <- apply(counts, 1, max) >= 10.0
keepCpm <- rowSums(cpm(counts)> 0.0) >= 0                  # Keep only genes with cpm above x in at least y samples

keep <- keepMin & keepCpm
x <- x[keep,]
counts <- counts[keep,]

# normalisation / voom
nf <- calcNormFactors(counts)
y<-voom(counts, design, plot=FALSE,lib.size=colSums(counts)*nf)

# prepare data for violin plot
# genelist<-read.csv('geneids.txt',header=F)
z0<-counts
z1<-log(cpm(counts)+1)
z2<-y$E
# rownames(z0)<-genelist$V1
# rownames(z1)<-genelist$V1
# rownames(z2)<-genelist$V1
rownames(z0)<-x$Geneid
rownames(z1)<-x$Geneid
rownames(z2)<-x$Geneid
write.table(z0,'raw_expression.tsv',sep="\t")
write.table(z1,'logcpm_expression.tsv',sep="\t")
write.table(z2,'normalised_expression.tsv',sep="\t")

# gene function description
d<-read.table('data_ab5075/GCF_000963815.1_ASM96381v1_feature_table.txt',sep="\t",quote="")
d1<-data.frame(locus=d[which(d$V1=='CDS'),c('V17')],gid=d[which(d$V1=='CDS'),c('V15')],product=d[which(d$V1=='CDS'),c('V15')])

# linear model fitting
fit <- lmFit(y,design)

out2<-data.frame()
degs<-data.frame()
for (i in 1:length(treatment.list)){
  cont.matrix1 <- subset(cont.matrix,select=i)
  fit2 <- contrasts.fit(fit, cont.matrix1)
  fit2 <- eBayes(fit2)
  out <- topTable(fit2, n=Inf, sort.by='none')
  colnames(out)<-paste(treatment.list[i],colnames(out),sep=".")
  if (i == 1){
    out2<-cbind(out)
    degs<-cbind(as.numeric(abs(out[,paste(treatment.list[i],'logFC',sep=".")])>=1 & out[,paste(treatment.list[i],'adj.P.Val',sep=".")]<=0.05))
  }else{
    out2<-cbind(out2,out)
    degs<-cbind(degs,as.numeric(abs(out[,paste(treatment.list[i],'logFC',sep=".")])>=1 & out[,paste(treatment.list[i],'adj.P.Val',sep=".")]<=0.05))
  }
}
out2<-cbind(out2,x[, export_cols])
# out2<-merge(out2,d1,by.x="Geneid",by.y="locus",all.x=T) # keeping all degs incl. RNA encoding genes
write.csv(out2,file='results_ab.csv')

colnames(degs)<-treatment.list
rownames(degs)<-out2$Geneid
# upset plot
degs<-data.frame(degs)

# tiff('ab_upset.low.1h.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
# upset(degs[,c(1,3,5,7)])
# dev.off()
# tiff('ab_upset.high.1h.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
# upset(degs[,c(2,4,6,8)])
# dev.off()
# tiff('ab_upset.lowhigh.1h.tiff', units="in", width=15, height=8, res=300, compression = 'lzw')
# upset(degs[,c(1:8)],nsets=8,nintersects=NA)
# dev.off()
#
# tiff('ab_upset.low.4h.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
# upset(degs[,c(9,11,13,15)])
# dev.off()
# tiff('ab_upset.high.4h.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
# upset(degs[,c(10,12,14,16)])
# dev.off()
# tiff('ab_upset.lowhigh.4h.tiff', units="in", width=15, height=8, res=300, compression = 'lzw')
# upset(degs[,c(9:16)],nsets=8,nintersects=NA)
# dev.off()
#
# tiff('ab_upset.all.tiff', units="in", width=30, height=8, res=300, compression = 'lzw')
# upset(degs[,c(1:16)],nsets=16,nintersects=NA)
# dev.off()
pdf('allsets.upset.pdf',height=30,width = 10)
degs.all0<-make_comb_mat(degs[rowSums(degs[,1:16])>=1,1:16], mode='distinct')
degs.all<-degs.all0[comb_size(degs.all0) > 1]
getPalette<-colorRampPalette(brewer.pal(11,'Spectral'))
upset.all<-UpSet(t(degs.all), pt_size=unit(1.5,'mm'),
           lwd=0.5,set_order=order(set_name(degs.all)),
           comb_col = getPalette(max(comb_degree(degs.all)))[comb_degree(degs.all)],
           top_annotation = columnAnnotation(
             "No. of DEGs" = anno_barplot(set_size(degs.all),
                                          axis_param = list(direction = "normal"),
                                          border = FALSE,
                                          gp = gpar(fill = "black"),
                                          width = unit(2, "cm"),height=unit(3,'cm'),cex=0.5
             )),
           row_names_side = "top",height = unit(30, "cm"),
           column_labels = paste0(rep(c('Colistin','F287','F319','F365'),each=2,2),'-',rep(c('L','H'),8),'-',c(rep('1h',8),rep('4h',8))),
           column_names_gp = gpar(fontsize = 6),width = unit(4, "cm"))
print(upset.all)
dev.off()
pdf('sets.L.upset.pdf',height=15,width = 10)
degs.0.L<-make_comb_mat(degs[rowSums(degs[,seq(1,16,2)])>=1,seq(1,16,2)], mode='distinct')
degs.L<-degs.0.L[comb_size(degs.0.L) > 1]
getPalette<-colorRampPalette(brewer.pal(11,'Spectral'))
upset.L<-UpSet(t(degs.L), pt_size=unit(3,'mm'),
                 lwd=1.5,set_order=1:8,
                 comb_col = getPalette(max(comb_degree(degs.L)))[comb_degree(degs.L)],
                 top_annotation = columnAnnotation(
                   "No. of DEGs" = anno_barplot(set_size(degs.L),
                                                axis_param = list(direction = "normal"),
                                                border = FALSE,
                                                gp = gpar(fill = "black"),
                                                width = unit(2, "cm"),height=unit(1,'cm'),cex=0.5
                   )),
                 row_names_side = "top",height = unit(30, "cm"),
                 column_labels = paste0(rep(c('Colistin','F287','F319','F365'),each=1,2),'-',rep('L',8),'-',c(rep('1h',4),rep('4h',4))),
                 column_names_gp = gpar(fontsize = 15),width = unit(6, "cm"),
                 right_annotation = upset_right_annotation(t(degs.L),width=unit(1,"cm"),
                                                           gp = gpar(fill = getPalette(max(comb_degree(degs.L)))[comb_degree(degs.L)]))
                )
print(upset.L)
dev.off()
pdf('sets.H.upset.pdf',height=20,width = 10)
degs.0.H<-make_comb_mat(degs[rowSums(degs[,seq(2,16,2)])>=1,seq(2,16,2)], mode='distinct')
degs.H<-degs.0.H[comb_size(degs.0.H) > 1]
getPalette<-colorRampPalette(brewer.pal(11,'Spectral'))
upset.H<-UpSet(t(degs.H), pt_size=unit(2,'mm'),
               lwd=1,set_order=1:8,
               comb_col = getPalette(max(comb_degree(degs.H)))[comb_degree(degs.H)],
               top_annotation = columnAnnotation(
                 "No. of DEGs" = anno_barplot(set_size(degs.H),
                                              axis_param = list(direction = "normal"),
                                              border = FALSE,
                                              gp = gpar(fill = "black"),
                                              width = unit(2, "cm"),height=unit(max(set_size(degs.H))/max(set_size(degs.L)),'cm'),cex=0.5
                 )),
               row_names_side = "top",height = unit(30, "cm"),
               column_labels = paste0(rep(c('Colistin','F287','F319','F365'),each=1,2),'-',rep('H',8),'-',c(rep('1h',4),rep('4h',4))),
               column_names_gp = gpar(fontsize = 15),width = unit(6, "cm"),
               right_annotation = upset_right_annotation(t(degs.H),width =unit(max(comb_size(degs.H))/max(comb_size(degs.L)),"cm"),
                                                         gp = gpar(fill = getPalette(max(comb_degree(degs.H)))[comb_degree(degs.H)]))
)
print(upset.H)
dev.off()
k<-apply(degs,2,sum)
write.csv(k,file='num_degs_ab.csv')
save(degs.all,degs.all0,file='upset_comb.RData')
save(degs.L,degs.0.L,file='upset_comb_L.RData')
#
#
# # Volcano plot
# volcanoplot(fit2,coef=1,names = x[,'Geneid'],highlight=sum(abs(dt[,'Col.L_Cont_1h'])))

# # venn chart
# # 1h
# # colistin L&H
# vennDiagram(dt[,1:2],circle.col=cols[1:2])
# # colistin L&H
# vennDiagram(dt[,3:4])
# # colistin L&H
# vennDiagram(dt[,5:6])
# # colistin L&H
# vennDiagram(dt[,7:8])


# tune pca components
pdf('tune.pca.ab.pdf')
tune.pca(t(log(cpm(counts)+1)))
dev.off()

# pca

s$Treatment <- paste(s$Antibiotic, s$Concentration, sep="_")
results.pca<-pca(t(log(cpm(counts)+1)),ncomp=2)
pdf('pca.ab.pdf')
plotIndiv(results.pca,group=s$Treatment,style='ggplot2',
  pch=as.factor(s$Time), legend=T,col=cols,alpha=0.2,
  legend.title='Treatment', legend.title.pch='Time (h)',ellipse=F)
dev.off()

# cal FPKM
y1<-DGEList(counts=data.matrix(counts),genes = x[,c('Geneid','Length')])
y1<-calcNormFactors(y1)
y_rpkm<-rpkm(y1,gene.length = as.numeric(y1$genes$Length))
rownames(y_rpkm)<-rownames(degs)
y_rpkm<-data.frame(y_rpkm)
save(y_rpkm,file='5075.rpkm.RData')

# splsda
if (0){
for (i in c(1,4)){
  for (j in c(2,8)){
    class<-as.factor(s[which(s$Time==i & s$Concentration==j ),'Antibiotic'])
    sample.id<-s[which(s$Time==i & s$Concentration==j ),'Sample.ID']
    gene.exp<-t(log(cpm(counts[,colnames(counts) %in% sample.id]+1)))

    results.splsda<-splsda(gene.exp,class,ncomp=3,keepX=c(100,100,100))
    tiff(paste(c('splsda.ab',i,j,'tiff'),sep="",collapse="."), units="in", width=15, height=8, res=300, compression = 'lzw')
    plotIndiv(results.splsda, ind.names=F,legend=T,ellipes=T,title='AB5075 RNASeq Col F287 F319 F365')
    dev.off()
    selectVar(results.splsda,comp=1)$value
    # plotLoadings(results.splsda,contrib='max',method='mean',comp=2,name.var=genelist$V1)
    tiff(paste(c('splsda.ab.loading','comp1',i,j,'tiff'),sep="",collapse="."), units="in", width=15, height=8, res=300, compression = 'lzw')
    plotLoadings(results.splsda,contrib='max',method='mean',comp=2,name.var=x$Geneid)
    dev.off()
  }
}
}


# check expression of interested gene
if (0){
  gene.name<-'ABUW_RS04050'
  expr<-data.frame(exp=as.numeric(z1[gene.name,]),
      treatment=as.factor(paste(s$Antibiotic, s$Concentration, s$Time,sep="_")))

  p<-ggplot(expr,aes(x=treatment,y=exp,color=treatment))+
      geom_jitter(position=position_jitter(0.2))+
      labs(title=paste0(gene.name,'_all'),x="Treatment", y = "Expression (log(cpm))")
  p + theme_classic()

  expr1<-data.frame(exp=as.numeric(z1[gene.name,colnames(z1) %in% sample.id]),
      treatment=class)

  p1<-ggplot(expr,aes(x=treatment,y=exp,color=treatment))+
      geom_jitter(position=position_jitter(0.2))+
      labs(title=gene.name,x="Treatment", y = "Expression (log(cpm))")
  p1 + theme_classic()
  ggsave('ABUW_RS04050.png')
}

# heatmap
# tiff("heatmap.norm.ab.tiff", units="in", width=8, height=5, res=300)
# pheatmap(mat = as.matrix(y$E),color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255),scale = "column")
#          # annotation_col = annotation_col,
#          # annotation_colors = ann_colors,
# # )
# dev.off()

save(y,degs,design,s,k,g,out2,treatment.list,file='degs.5075.Rdata')
