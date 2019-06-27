
cat ../annot/annot.gtf | cut -f1,4,5 | sort | uniq > annot.3

soft_fasta_to_3 < ../data/genome.fasta.masked > mask.3
cat mask.3 | awk '{if($3-$2+1>=50) print }' > mask_50.3
cat mask.3 | awk '{if($3-$2+1>=100) print }' > mask_100.3
cat mask.3 | awk '{if($3-$2+1>=500) print }' > mask_500.3
cat mask.3 | awk '{if($3-$2+1>=1000) print }' > mask_1000.3

bedtools intersect  -a annot.3 -b mask.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b mask_50.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b mask_100.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b mask_500.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b mask_1000.3 -wo | cut -f1,2,3 | sort | uniq | wc -l

cat mask.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat mask_50.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat mask_100.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat mask_500.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat mask_1000.3 | awk '{sum+=$3-$2+1}END{print sum}'

