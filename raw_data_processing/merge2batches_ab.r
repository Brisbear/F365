library(ggplot2)
library(tidyverse)

df<-read.table('combinedCounts_ab_final.txt',sep="\t",header=T,quote="")
tobekept<-unique(gsub(".*raw_data.(.*)_L.*",'\\1',colnames(df),perl=T))
df1<-df[,c(1:62,64)]
df1[,62]<-df[,62]+df[,63]
colnames(df1)<-tobekept

write.table(df1,'merged.ab.counts.txt',sep="\t",quote=F)
