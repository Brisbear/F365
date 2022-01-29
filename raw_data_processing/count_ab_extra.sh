#! /bin/bash

# featureCounts -p -T 30 -t gene -g gene_id -a Ab_AB5075.gtf -o D-4h-8_L1_1.pe.extra.bam.counts.txt ../N2003266_JXZ_30-372854748/raw_data/D-4h-8_L1_1.pe.bam  
featureCounts -p -T 30 -t gene -g gene_id -a GCF_000963815.1_ASM96381v1_genomic.gtf -o D-4h-8_L1_1.pe.extra.bam.counts.txt ../N2003266_JXZ_30-372854748/raw_data/D-4h-8_L1_1.pe.bam  


