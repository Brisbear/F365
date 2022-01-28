library(ComplexHeatmap)
library(RColorBrewer)

load('degs.hk2.Rdata')

a<-out2[,grep('logFC',colnames(out2))]
up<- 1* as.data.frame(a>0 & degs)
down<- 1* as.data.frame(a<0 & degs)


degs4=make_comb_mat(PMB.up=out2$Geneid[up$PMB==1],
                    PMB.down=out2$Geneid[down$PMB==1],
                    Colistin.up=out2$Geneid[up$Colistin==1],
                    Colistin.down=out2$Geneid[down$Colistin==1],
                    F319.up=out2$Geneid[up$F319==1],
                    F319.down=out2$Geneid[down$F319==1],
                    F287.up=out2$Geneid[up$F287==1],
                    F287.down=out2$Geneid[down$F287==1],
                    F365.up=out2$Geneid[up$F365==1],
                    F365.down=out2$Geneid[down$F365==1],
                    mode='distinct')
degs5=make_comb_mat(PMB.up=out2$Geneid[up$PMB==1],
                    PMB.down=out2$Geneid[down$PMB==1],
                    Colistin.up=out2$Geneid[up$Colistin==1],
                    Colistin.down=out2$Geneid[down$Colistin==1],
                    F319.up=out2$Geneid[up$F319==1],
                    F319.down=out2$Geneid[down$F319==1],
                    F287.up=out2$Geneid[up$F287==1],
                    F287.down=out2$Geneid[down$F287==1],
                    F365.up=out2$Geneid[up$F365==1],
                    F365.down=out2$Geneid[down$F365==1],
                    mode='intersect')

pdf('ht.up_down.upset.distinct.pdf',width = 12)
temp.col<-brewer.pal(12,"Paired")
ht1<-UpSet(degs4, pt_size=unit(3,'mm'),
           lwd=1,comb_col = brewer.pal(8,"Set2")[comb_degree(degs4)],
           set_order = 1:10,
           right_annotation = rowAnnotation(
             "No. of DEGs" = anno_barplot(set_size(degs4),
                                          axis_param = list(direction = "normal"),
                                          border = FALSE,
                                          gp = gpar(fill = rep(c('grey','black'),5)),
                                          width = unit(2, "cm")
             )),row_names_side = "left",height = unit(5, "cm"),width = unit(15, "cm"))


ht2<-UpSet(degs5, pt_size=unit(3,'mm'),
           lwd=1,comb_col = brewer.pal(8,"Set2")[comb_degree(degs5)],
           set_order = 1:10,
           right_annotation = rowAnnotation(
             "No. of DEGs" = anno_barplot(set_size(degs5),
                                          axis_param = list(direction = "normal"),
                                          border = FALSE,
                                          gp = gpar(fill = rep(c('grey','black'),5)),
                                          width = unit(2, "cm")
             )),row_names_side = "left",height = unit(5, "cm"),width = unit(15, "cm"))
print(ht1)
# print(ht2)
dev.off()
