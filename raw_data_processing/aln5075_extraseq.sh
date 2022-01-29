#! /bin/bash

subread-align -t 0 -T 30 -d 50 -D 600 -i ../00_ref/Ab_AB5075  -r ../N2003266_JXZ_30-372854748/raw_data/D-4h-8_L1_1.fq.gz -R ../N2003266_JXZ_30-372854748/raw_data/D-4h-8_L1_2.fq.gz -o ../N2003266_JXZ_30-372854748/raw_data/D-4h-8_L1_1.pe.bam --quality 2>&1 | tee ../N2003266_JXZ_30-372854748/raw_data/D-4h-8_L1_1.pe.log


