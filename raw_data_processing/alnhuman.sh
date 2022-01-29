#! /bin/bash

for f in ../raw_data/A*1.fq.gz; do
	echo $f
	base=${f%.fq.gz}
	log=$base".pe.log"
	output=$base".pe.bam"
	F=${f/1.fq.gz/2.fq.gz}
	echo $F
#	subread-align -t 0 -T 30 -d 50 -D 600 -i ../00_ref/Ab_AB5075  -r $f -R $F -o $output --quality 2>&1 | tee $log
	subjunc  -T 30 -d 50 -D 600 -i ../00_ref/GRCh38 -r $f -R $F -o $output --quality 2>&1 | tee $log
done




