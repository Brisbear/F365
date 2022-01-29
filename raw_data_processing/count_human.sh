#! /bin/bash


for f in ../raw_data/A*.bam; do
	base=${f##.*\/}
	output=$base".counts.txt"
	featureCounts -T 30 -t exon -g gene_id -a Homo_sapiens.GRCh38.94.gtf -o $output $f
done




 
