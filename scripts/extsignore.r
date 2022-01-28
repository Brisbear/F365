library(ggplot2)
library(igraph)
library(RColorBrewer)
library(tidyverse)
# library(UpSetR)
library(ComplexHeatmap)
gen.subgraph<-function(net.all,genes,condition,abslogfc,color,t){
  # print(length(which(V(net.all)$name %in% genes)))
  # print(genes)
  V(net.all)$color[which(V(net.all)$name %in% genes)]<-color
  V(net.all)$size[which(V(net.all)$name %in% genes)]<-abslogfc
  selnodes <- V(net.all)[name %in% genes]
  selegoV <- ego(net.all, order=1, nodes=selnodes,mode='all')
  selegoG <- induced_subgraph(net.all,unlist(selegoV))
  # print(V(selegoG)$name)
  # print(V(selegoG)$color)
  pdf(paste0('signor.net.',condition,'.pdf'),height=20,width = 20)
  # print(length(V(selegoG)$name))
  l1<-layout.auto(selegoG)
  l1 <- norm_coords(l1, ymin=-1, ymax=1, xmin=-1, xmax=1)
  # plot(selegoG,layout=l1)
  plot(selegoG,layout=l1)
  legend("bottomleft",legend = t$type,pch=21, col=t$color, 
         pt.bg=t$color,pt.cex=3, cex=2, ncol=1)
  dev.off()
  return(1)
}

load('degs.hk2.Rdata')
load('data_hk2/gene.alias.RData')
load('data_hk2/distinctive.upset.degs.RData')
load('data_hk2/degs.upset.type.RData')
network0<-read.table('data_hk2/signor_all_data_05_07_20.tsv',sep="\t",header=T,quote="",fill=T)
network <- network0 %>% filter(str_detect(MECHANISM,'transcriptional')) %>% 
  filter(!grepl('post',MECHANISM)) %>%
  filter(TYPEA=='protein' & TYPEB=='protein')
nodes<-data.frame(item=c(network$ENTITYA,network$ENTITYB),type=c(network$TYPEA,network$TYPEB)) %>%
  distinct(item,type) %>%
  filter(item!='') %>% 
  filter(!(item=='SRC' & type=='proteinfamily'))
types<-sort(unique(c(network$TYPEA,network$TYPEB)))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

types<-rbind(data.frame(type=types,color=getPalette(length(types))),
             data.frame(type=c('upregulated.genes','downregulated.genes'),color=c('red2','green2')))
types[types$type=='protein','color']<-'black'
b<-degs %>% filter(rowSums(degs[,1:5])>=1)
b$gene<-rownames(b)
out2[,grep('logFC',colnames(out2))][degs==0]<-0
e<-out2 %>% filter(rowSums(degs[,1:5])>=1)%>% 
  dplyr::select(contains('logFC'))
colnames(e)<-gsub('.logFC','',colnames(e))
e$gene<-rownames(e)

v.all<-data.frame(node=unique(c(network$ENTITYA,network$ENTITYB)))
network$DIRECT<-recode(network$DIRECT,UNK='NO')
e.all<-network %>% group_by(ENTITYA,ENTITYB,TYPEA,TYPEB,MECHANISM,EFFECT,DIRECT) %>%
  summarise(SCORE=sum(SCORE))
colnames(e.all)[1:2]<-c('from','to')
net.all<-graph_from_data_frame(d=e.all,vertices = v.all,directed=F)
node.color<-left_join(data.frame(item=V(net.all)$name),nodes,by='item') %>%
  left_join(types,by='type') %>% left_join(g.upset,by=c('item'='Gene.Symbol'))
node.color$col[is.na(node.color$col)]<-'black'
V(net.all)$color<-node.color$col
V(net.all)$frame.color<-NA
V(net.all)$size<-1
V(net.all)$label.color <- 'black'
V(net.all)$label.dist <- 0.3
V(net.all)$label.family	<-'Helvetica'
V(net.all)$label.cex <- 0.6
E(net.all)$width<-E(net.all)$SCORE
# E(net.all)$width<-0.5
E(net.all)$color[grepl('up',E(net.all)$EFFECT)]<-'#FF7F50'
E(net.all)$color[grepl('down',E(net.all)$EFFECT)]<-'#008B8B'
E(net.all)$color[!grepl('(up)|(down)',E(net.all)$EFFECT)]<-'#808080'
E(net.all)$arrow.size <- 0.4
E(net.all)$arrow.width <- 0.4
E(net.all)$arrow.mode<-E(net.all)$DIRECT %>% recode(YES=2,NO=0)
E(net.all)$curved <- 0.1


net.all.1<-net.all
V(net.all.1)$label<-NA
# pdf('signor.network.pdf',width=35,height=35)
# l <- layout_in_circle(net.all.1)
# l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
# plot(net.all.1,rescale=F, layout=l*1)
# dev.off()

# degs for any condition, ensembl gene column, filter out rows without gene symbol (noncoding genes and genes without assigned symbol)
v.sub0 <- e %>% left_join(gene.alias,by=c('gene'='Ensembl.Gene.ID')) %>% filter(!is.na(Gene.Symbol))
v.sub1<-intersect(v.sub0$Gene.Symbol, V(net.all)$name)

colref.merged.conditions<-unique(node.color[which(node.color$item %in% v.sub1),c('col','group')])%>% 
  arrange(match(group,comb_name(degs1)))
colnames(colref.merged.conditions)<-c('color','type')
gen.subgraph(net.all,v.sub1,'anytreatment',rep(4,length(v.sub1)),
             node.color[which(node.color$item %in% v.sub1),'col'],colref.merged.conditions) 
# subgraph by shortest paths between each possible pair of nodes
# sp<-all_shortest_paths(net.all,from = v.sub1,to = v.sub1,mode ='all')
# v.sub2<-sort(unique(unlist(sp$res)))
# net.sub<-induced_subgraph(net.all,V(net.all)$name[v.sub2])
# plot(net.sub,layout=layout_with_graphopt(net.sub))


for (i in 1:5){
  genes<-e %>% filter(!!sym(colnames(b)[i])!=0) %>%
    dplyr::select(!!sym(colnames(b)[i]),gene) %>%
    left_join(gene.alias,by=c('gene'='Ensembl.Gene.ID')) %>%
    filter(!is.na(Gene.Symbol)) %>%
    filter(Gene.Symbol %in% V(net.all)$name)
  print(length(genes$Gene.Symbol))
  c<-rep('red2',nrow(genes))
  c[genes[,1]<0]<-'green2'
  s<-abs(genes[,1])*1.5
  gen.subgraph(net.all,genes$Gene.Symbol,colnames(b)[i],s,c,types)
}

