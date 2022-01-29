#! /bin/bash

for f in ../raw_data/*-*.bam; do
	base=${f##.*\/}
	output=$base".counts.txt"
	featureCounts -p -T 30 -t gene -g gene_id -a Ab_AB5075.gtf -o $output $f

done
 
