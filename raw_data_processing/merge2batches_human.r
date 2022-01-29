library(ggplot2)
library(tidyverse)

df<-read.table('combinedCounts_human_final.txt',sep="\t",header=T,quote="")
tobekept<-paste0('A',c(3:6,8,10,13:17,19,23:27,30))
df1<-cbind(df[,c(1:6)], select(df,contains(tobekept)))
df2<-df[,c(1:6)]
for (i in 1:length(tobekept)){
	a<-select(df1,contains(paste0(tobekept[i],'_')))
	b<-data.frame(x=rowSums(a))
	colnames(b)<-tobekept[i]
	df2<-cbind(df2,b)
}

write.table(df2,'merged.hk2.counts.txt',sep="\t",quote=F)
