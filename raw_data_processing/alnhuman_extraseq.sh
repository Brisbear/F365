#! /bin/bash

for f in ../N2003266_JXZ_30-372854748/raw_data/A*1.fq.gz; do
	echo $f
	base=${f%.fq.gz}
	log=$base".pe.log"
	output=$base".pe.bam"
	base=${base%1}
	F=$base"2.fq.gz"
	echo $F
	subjunc  -T 30 -d 50 -D 600 -i ../00_ref/GRCh38 -r $f -R $F -o $output --quality 2>&1 | tee $log
done




