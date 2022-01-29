library(tidyverse)

df<-read.csv('results_ab.csv')

exp<- df %>% select(matches("Geneid") | matches(".*(logFC)|(adj.P)"))

comparisons<-colnames(df)[grepl('logFC',colnames(df))]
comparisons<-lapply(comparisons,function(x) {gsub('.logFC','',x)})

if (1){
  for (i in 1:length(comparisons)){
    comp <- exp %>% select(matches("Geneid") | contains(comparisons[[i]]))
    # degs
    comp.de <- comp[which(abs(comp[,2])>=1 & comp[,3] <=0.05),]
    write.table(comp.de,paste(c(comparisons[[i]],'biocyc','txt'),sep="",collapse="."),row.names=F,quote=F,col.names=F,sep="\t")
  }
}


keep<- apply(abs(exp %>% select(contains('logFC'))),1,function(x) any(x>=1)) &
    apply(abs(exp %>% select(contains('adj.P'))),1,function(x) any(x<=0.05))

exp.fc<-exp[keep,] %>% select(contains('logFC'))
exp.p<-exp[keep,] %>% select(contains('adj.P'))

exp.fc[abs(exp.fc)<1 | exp.p>0.05]<- 0

a<-exp[keep,1]
exp.fc<-cbind(a,exp.fc)
colnames(exp.fc)[1]<-'Geneid'

write.table(exp.fc,'combined_logFC_ab.txt',sep="\t",row.names=F,quote=F)
