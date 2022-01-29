#! /bin/bash


for f in ../N2003266_JXZ_30-372854748/raw_data/A*.bam; do
	base=${f##.*\/}
	output=$base".extra.counts.txt"
	featureCounts -T 30 -t exon -g gene_id -a Homo_sapiens.GRCh38.94.gtf -o $output $f
done




 
