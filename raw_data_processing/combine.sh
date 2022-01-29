#! /bin/bash


# The output file
#out="combinedCounts_ab.txt"
out="combinedCounts_ab_final.txt"
# out="combinedCounts_human.txt"
# out="combinedCounts_human_final.txt"

# Combine all counts together into one doc
paste *-*.txt | tail -n +2 > tmpResults.txt
#paste A*.txt | tail -n +2 > tmpResults.txt

# Total number of columns
ncol=`awk '{print NF}' tmpResults.txt | sort -nu | tail -n 1`

# Extract the counts
awk -v f=1 -v t=$ncol -v m=7 -F "\t" '{	for(i=f; i<t ;i++) printf("%s",(i<m || i%m == 0 ? $i"\t" : "")); printf("%s\n", $t);}' tmpResults.txt > $out

sed "s/\=/_/g" $out > tmpResults.txt
sed "1s/\-/_/g" tmpResults.txt > $out

# Remove the tmp file
rm tmpResults.txt

