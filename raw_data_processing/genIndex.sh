#! /bin/bash

# this uses old annotation of AB5075 genome
# subread-buildindex -o "Ab_AB5075" -B Ab_AB5075.fna

# this is new
subread-buildindex -o "Ab_AB5075" -B GCF_000963815.1_ASM96381v1_genomic.fna

# human genome reference was cipied from Harry PAO1/A549 infection project
