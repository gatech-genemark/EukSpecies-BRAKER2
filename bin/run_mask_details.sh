# cat genome.fasta.out | sed 's/NC_004354.4/X/' | sed 's/NT_033779.5/2L/' | sed 's/NT_033778.4/2R/' | sed 's/NT_037436.4/3L/' | sed 's/NT_033777.3/3R/' | sed 's/NC_004353.4/4/' | sed 's/NC_024512.1/Y/' > z
# mv z genome.fasta.out 

# cat genome.fasta.out | sed 's/NC_003070.9/Chr1/' | sed 's/NC_003070.9/Chr1/' | sed 's/NC_003074.8/Chr3/' | sed 's/NC_003075.7/Chr4/' | sed 's/NC_003076.8/Chr5/' > z
# mv z genome.fasta.out 

cat genome.fasta.out | grep '_family-' | grep ' Unknown '    | awk '{print $5 "\t" $6 "\t" $7}' > Unknown.3
cat genome.fasta.out | grep '_family-' | grep -v ' Unknown ' | awk '{print $5 "\t" $6 "\t" $7}' > Known.3
tail -n +3 genome.fasta.out | grep -v '_family-'             | awk '{print $5 "\t" $6 "\t" $7}' > Low.3

cat Low.3 | awk '{if($3-$2+1>=50) print }'   > Low_50.3
cat Low.3 | awk '{if($3-$2+1>=100) print }'  > Low_100.3
cat Low.3 | awk '{if($3-$2+1>=500) print }'  > Low_500.3
cat Low.3 | awk '{if($3-$2+1>=1000) print }' > Low_1000.3

cat Unknown.3 | awk '{if($3-$2+1>=50) print }'   > Unknown_50.3
cat Unknown.3 | awk '{if($3-$2+1>=100) print }'  > Unknown_100.3
cat Unknown.3 | awk '{if($3-$2+1>=500) print }'  > Unknown_500.3
cat Unknown.3 | awk '{if($3-$2+1>=1000) print }' > Unknown_1000.3

cat Known.3 | awk '{if($3-$2+1>=50) print }'   > Known_50.3
cat Known.3 | awk '{if($3-$2+1>=100) print }'  > Known_100.3
cat Known.3 | awk '{if($3-$2+1>=500) print }'  > Known_500.3
cat Known.3 | awk '{if($3-$2+1>=1000) print }' > Known_1000.3

echo "# num overlap"
echo "Low"
bedtools intersect  -a annot.3 -b Low.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Low_50.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Low_100.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Low_500.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Low_1000.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
echo "Unknown"
bedtools intersect  -a annot.3 -b Unknown.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Unknown_50.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Unknown_100.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Unknown_500.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Unknown_1000.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
echo "Known"
bedtools intersect  -a annot.3 -b Known.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Known_50.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Known_100.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Known_500.3 -wo | cut -f1,2,3 | sort | uniq | wc -l
bedtools intersect  -a annot.3 -b Known_1000.3 -wo | cut -f1,2,3 | sort | uniq | wc -l

echo "# sum overlap"
echo "Low"
cat Low.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Low_50.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Low_100.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Low_500.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Low_1000.3 | awk '{sum+=$3-$2+1}END{print sum}'
echo "Unknown"
cat Unknown.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Unknown_50.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Unknown_100.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Unknown_500.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Unknown_1000.3 | awk '{sum+=$3-$2+1}END{print sum}'
echo "Known"
cat Known.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Known_50.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Known_100.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Known_500.3 | awk '{sum+=$3-$2+1}END{print sum}'
cat Known_1000.3 | awk '{sum+=$3-$2+1}END{print sum}'

echo "# count"
echo "Low"
wc -l Low.3
wc -l Low_50.3
wc -l Low_100.3
wc -l Low_500.3
wc -l Low_1000.3
echo "Unknown"
wc -l Unknown.3
wc -l Unknown_50.3
wc -l Unknown_100.3
wc -l Unknown_500.3
wc -l Unknown_1000.3
echo "Known"
wc -l Known.3
wc -l Known_50.3
wc -l Known_100.3
wc -l Known_500.3
wc -l Known_1000.3


