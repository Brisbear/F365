library(limma)
library(edgeR)
library(jsonlite)
library(mixOmics)
library(RColorBrewer)
library(ComplexHeatmap)
library(UpSetR)

load("data_hk2/GRCh38.94.gtf.CDS.simple.RData")
cols <- brewer.pal(6, "Dark2")
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
counts_file <- 'data_hk2/merged.hk2.counts.txt'
sample.list<-read.table('data_hk2/sample.final.tsv',sep="\t",header=T,quote="")
count_cols <-sample.list$Sample.ID

design <- matrix(c(as.numeric(grepl('Control',sample.list$Treatment)),
                   as.numeric(grepl('PMB',sample.list$Treatment)),
                   as.numeric(grepl('Colistin',sample.list$Treatment)),
                   as.numeric(grepl('F319',sample.list$Treatment)),
                   as.numeric(grepl('F287',sample.list$Treatment)),
                   as.numeric(grepl('F365',sample.list$Treatment))), ncol=6,
dimnames=list(count_cols,c('Control','PMB','Colistin','F319','F287','F365')))


cont.matrix <- matrix(c(c(-1,1,0,0,0,0),c(-1,0,1,0,0,0),c(-1,0,0,1,0,0),c(-1,0,0,0,1,0),c(-1,0,0,0,0,1)), ncol=5, dimnames=list(c('Control','PMB','Colistin','F319','F287','F365'),c('PMB','Colistin','F319','F287','F365')))

export_cols <- c(c('Geneid','Chr','Start','End','Strand','Length'),
              count_cols[c(grep('Control',sample.list$Treatment),grep('PMB',sample.list$Treatment),
              grep('Colistin',sample.list$Treatment),grep('F319',sample.list$Treatment),
              grep('F287',sample.list$Treatment),grep('F365',sample.list$Treatment))])

exp_cols_network <- c('Geneid',count_cols[c(grep('Control',sample.list$Treatment),grep('PMB',sample.list$Treatment),
grep('Colistin',sample.list$Treatment),grep('F319',sample.list$Treatment),
grep('F287',sample.list$Treatment),grep('F365',sample.list$Treatment))])

x<-read.delim(counts_file,
              sep="\t",
              check.names=FALSE,
              colClasses='character',
              na.strings=c(),
              skip=0
              )
rownames(x)<-x$Geneid
# Now re-read the first header line.  Workaround R problem that always has strip.white=T for read.table
colnames(x) <- scan(counts_file,
                    what="",
                    sep="\t",
                    nlines=1,
                    strip.white=F,
                    quote = "\"",
                    skip="0"
                    )

x[,count_cols] <- apply(x[,count_cols], 2, function(v) as.numeric(v))     # Force numeric count columns
counts <- x[, count_cols]

keepMin <- apply(counts, 1, max) >= 10.0
keepCpm <- rowSums(cpm(counts)> 0.0) >= 0                  # Keep only genes with cpm above x in at least y samples

# keepCol<- !col(counts[1,]) %in% c(11:12, 14, 17, 21:27, 30)
keep <- keepMin & keepCpm

x <- x[keep,]
counts <- counts[keep,]
# x <- x[keep,c(rep(TRUE,6),keepCol)]
# counts <- counts[keep,keepCol]
# design1<-design[keepCol,]

# cal FPKM 
y1<-DGEList(counts=data.matrix(counts),genes = x[,c('Geneid','Length')])
y1<-calcNormFactors(y1)
y_rpkm<-rpkm(y1,gene.length = as.numeric(y1$genes$Length))
save(y_rpkm,file='hk2.rpkm.RData')


nf <- calcNormFactors(counts)
# y<-voom(counts, design1, plot=FALSE,lib.size=colSums(counts)*nf)
y<-voom(counts, design, plot=FALSE,lib.size=colSums(counts)*nf)


# prepare data for violin plot

genelist<-rownames(counts)
write.table(genelist,'geneids.txt',sep='\t',quote=F)
z0<-counts
z1<-log(cpm(counts)+1)
z2<-y$E
rownames(z0)<-genelist
rownames(z1)<-genelist
rownames(z2)<-genelist
write.table(z0,'raw_expression.tsv',sep="\t",quote=F)
write.table(z1,'logcpm_expression.tsv',sep="\t",quote=F)
write.table(z2,'normalised_expression.tsv',sep="\t",quote=F)


# fit <- lmFit(y,design1)
fit <- lmFit(y,design)

fc.threshold<-1.5

treatment.list<-c('PMB','Colistin','F319','F287','F365')
out2<-data.frame()
degs<-data.frame()
for (i in 1:5){
  cont.matrix1 <- subset(cont.matrix,select=i)
  fit2 <- contrasts.fit(fit, cont.matrix1)
  fit2 <- eBayes(fit2)
  out <- topTable(fit2, n=Inf, sort.by='none')
  colnames(out)<-paste(treatment.list[i],colnames(out),sep=".")
  if (i == 1){
    out2<-cbind(out)
    degs<-cbind(as.numeric(abs(out[,paste(treatment.list[i],'logFC',sep=".")])>=log(fc.threshold,2) & out[,paste(treatment.list[i],'adj.P.Val',sep=".")]<=0.05))
  }else{
    out2<-cbind(out2,out)
    degs<-cbind(degs,as.numeric(abs(out[,paste(treatment.list[i],'logFC',sep=".")])>=log(fc.threshold,2) & out[,paste(treatment.list[i],'adj.P.Val',sep=".")]<=0.05))
  }
}
out2<-cbind(out2,x[, export_cols])
write.table(out2,file='results_hk2.txt',sep="\t",quote=F)

out3<-x[,exp_cols_network]
colnames(out3)[1]<-'#NAME'

out3[seq(2,nrow(out3)+1),] <- out3[seq(1,nrow(out3)),]
out3[1,] <- c('#CLASS:group',sample.list[match(exp_cols_network[2:19], sample.list$Sample.ID),'Treatment'])
write.table(out3,file='results_hk2_for_network.txt',sep="\t",quote=F,row.names=F)

colnames(degs)<-treatment.list
rownames(degs)<-out2$Geneid
# upset plot
degs<-data.frame(degs)
# degs1=make_comb_mat(PMB=out2$Geneid[degs$PMB==1],Colistin=out2$Geneid[degs$Colistin==1],F319=out2$Geneid[degs$F319==1],
#       F287=out2$Geneid[degs$F287==1],F365=out2$Geneid[degs$F365==1], mode='distinct',complement_size = comb_size(degs3)[1])
degs1=make_comb_mat(PMB=out2$Geneid[degs$PMB==1],Colistin=out2$Geneid[degs$Colistin==1],
      F287=out2$Geneid[degs$F287==1],F319=out2$Geneid[degs$F319==1],F365=out2$Geneid[degs$F365==1], mode='distinct')
degs2=make_comb_mat(PMB=out2$Geneid[degs$PMB==1],Colistin=out2$Geneid[degs$Colistin==1],
      F287=out2$Geneid[degs$F287==1],F319=out2$Geneid[degs$F319==1],F365=out2$Geneid[degs$F365==1], mode='intersect')
degs3=make_comb_mat(PMB=out2$Geneid[degs$PMB==1],Colistin=out2$Geneid[degs$Colistin==1],
      F287=out2$Geneid[degs$F287==1],F319=out2$Geneid[degs$F319==1],F365=out2$Geneid[degs$F365==1], mode='union')

save(degs1,file='distinctive.upset.degs.RData')
# prepare heatmaps
pdf('ht1.upset.pdf')
ht1<-UpSet(degs1, pt_size=unit(5,'mm'),
    lwd=1,
    # comb_col = brewer.pal(9,"Set1")[comb_degree(degs1)+1],
    comb_col = getPalette(length(comb_size(degs1)))[1:length(comb_size(degs1))],
    right_annotation = rowAnnotation(
    "No. of DEGs" = anno_barplot(set_size(degs1),
        axis_param = list(direction = "normal"),
        border = FALSE,
        gp = gpar(fill = "black"),
        width = unit(2, "cm")
    )),row_names_side = "left",height = unit(2, "cm"),width = unit(10, "cm"))
print(ht1)
dev.off()
pdf('ht2.upset.pdf')
ht2<-UpSet(degs2, pt_size=unit(5,'mm'),
    lwd=1,comb_col = brewer.pal(9,"Set1")[comb_degree(degs2)+1],
    right_annotation = rowAnnotation(
    "No. of DEGs" = anno_barplot(set_size(degs2),
    axis_param = list(direction = "normal"),
    border = FALSE,
    gp = gpar(fill = "black"),
    width = unit(2, "cm")
  )),row_names_side = "left",height = unit(2, "cm"),width = unit(12, "cm"))
print(ht2)
dev.off()
# tiff('hk2_upset.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
tiff('hk2_upset.tiff', units="in", width=10, height=8, res=100, compression = 'lzw')
upset(degs,keep.order=T,order.by='degree',group.by = "degree",point.size = 5,sets.x.label = "No. of DEGs",text.scale = 1.5,
    queries = list(list(query = intersects,params = list('F365','PMB','Colistin','F287','F319'), color = "red", active = T),
              list(query = intersects,params = list('PMB'), color = "blue", active = T),
              list(query = intersects,params = list('F365'), color = "green", active = T),
              list(query = intersects,params = list('PMB','Colistin'), color = "orange", active = T)))

dev.off()
k<-data.frame(apply(degs,2,sum))
write.csv(k,file='num_degs_hk2.csv')

# pca

# s<-read.table('sample_reduced.tsv',sep="\t",header=T)

# results.pca<-pca(t(y$E),ncomp=10)
# e <- DGEList(counts=data.matrix(counts), genes=x[,c("Geneid","Length")])
# e_rpkm <- rpkm(e, gene.length=as.numeric(e$genes$Length))
# results.pca<-pca(t(log(e_rpkm+min(e_rpkm[e_rpkm!=0])+1)),ncomp=2)
results.pca<-pca(t(log(cpm(counts)+1)),ncomp=2)
pdf('pca_hk2.pdf')
# plotIndiv(results.pca,group=s$Treatment,style='ggplot2',legend=T,col=cols,pch=21,
#     legend.title='Treatment', title='HK2 RNASeq control - PMB - Colistin - F287 - F319 - F365')
plotIndiv(results.pca,group=sample.list$Treatment,style='graphics',
      legend=T,col=cols,pch=21,ellipse=T,legend.title='Treatment',title='PCA HK-2')
dev.off()

# tiff("heatmap.norm.hk2.tiff", units="in", width=8, height=5, res=300)
# pheatmap(mat = as.matrix(y$E),color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255),scale = "column")
#          # annotation_col = annotation_col,
#          # annotation_colors = ann_colors,
# dev.off()
degs.cds<-degs[which(rownames(degs) %in% cds$V9),]
out2.cds<-out2[which(rownames(out2) %in% cds$V9),]
save(degs,k,out2,degs.cds,out2.cds,treatment.list,file='degs.hk2.Rdata')

# plot the expressions using scattering points
b<-data.frame(logFC=numeric(),AveExpr=numeric(),treatment=character(),degs=character())
for (i in treatment.list){
  a<-out2[,paste0(i,c('.logFC','.AveExpr'))]
  colnames(a)<-c('logFC','AveExpr')
  a$treatment<-i
  a$degs<-as.character(degs[,i])
  b<-rbind(b,a)
}

p<-ggplot(b,aes(logFC,AveExpr,colour=treatment,fill=degs))+geom_point(alpha=0.2,size=0.5)+theme_bw()
ggsave('rnaseqScatterplotHK2.png')