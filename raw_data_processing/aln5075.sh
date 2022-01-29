#! /bin/bash

for f in ../raw_data/*h*1.fq.gz; do
	base=${f%.fq.gz}
	log=$base".pe.log"
	output=$base".pe.bam"
	base=${base%1}
	F=$base"2.fq.gz"
	subread-align -t 0 -T 30 -d 50 -D 600 -i ../00_ref/Ab_AB5075  -r $f -R $F -o $output --quality 2>&1 | tee $log
done




