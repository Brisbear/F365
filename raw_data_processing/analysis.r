library(Rsubread)
fq.f <- list.files(path = "../raw_data/", pattern = '(B|C|D).*_1.fq.gz$',full.names = T)
fq.r <- list.files(path = "../raw_data", pattern = '(B|C|D).*_2.fq.gz$',full.names = T)
# buildindex(basename="ab5075uw",reference="GCF_000963815.1_ASM96381v1_genomic.fna")

for ( i in 1:length(fq.f)){
  print(i)
  align(index="ab5075uw",
        readfile1=fq.f[i],
        readfile2=fq.r[i],
        type='rna',
        nthreads=24)
}
# 
bam.files <- list.files(path = "../raw_data", 
                        pattern = "(B|C|D).*BAM$", 
                        full.names = T)
# props <- propmapped(files=bam.files)
# qs <- qualityScores(filename="raw_data/SN-RNA001-LCM6195_L2_1.fq.gz",nreads=100)
# boxplot(qs)
fc <- featureCounts(bam.files, 
                    annot.ext = 'GCF_000963815.1_ASM96381v1_genomic1.gtf',
                    isGTFAnnotationFile = T,
                    GTF.featureType='tmRNA',
                    GTF.attrType = "gene_id",
                    nthreads = 24,
                    isPairedEnd=T,
                    genome='GCF_000963815.1_ASM96381v1_genomic.fna',
                    useMetaFeatures = F,
                    verbose = T)
save(fc,file='counts.Rdata')
