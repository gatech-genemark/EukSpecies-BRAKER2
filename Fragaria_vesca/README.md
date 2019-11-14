# Species: _Fragaria vesca  Assembly 4_  
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup  
```
species="Fragaria_vesca_4"
base="/storage3/EukSpecies"
export PATH="$base/bin:$PATH"
export base="$base/$species"
cd $base
if [ "$(pwd)" != "$base" ]; then echo "error, folder not found: $base"; fi
umask 002
```
Create core folders  
```
cd $base
mkdir -p arx annot data
```
### Genome sequence  
Assembly description is at 
```
# download data
cd $base/arx
wget ftp://ftp.bioinfo.wsu.edu/www.rosaceae.org/Fragaria_vesca/Fvesca-genome.v4.0.a1/assembly/Fragaria_vesca_v4.0.a1.fasta.gz
gunzip Fragaria_vesca_v4.0.a1.fasta.gz
flip -u Fragaria_vesca_v4.0.a1.fasta

# create ID table
grep '^>' F*.fasta > deflines
# move sequence IDs to "MtrunA17Chr" style
cat deflines | cut -b2- | awk '{print $1 "\t" $1}' > list.tbl

# select and reformat sequence; all uppercase 
get_fasta_with_tag.pl --swap --in F*.fasta  --out tmp_genome.fasta  --list list.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

# check sequence and clean folder
probuild --stat --details --seq genome.fasta
probuild --stat --details --seq tmp_genome.fasta

# put in work folder
mv genome.fasta ../data/genome.fasta

# clean tmp files
rm tmp_genome.fasta
gzip F*.fasta
```
### Masking: _de novo_ and _species specific_
Run _de novo_ masking of genome using RepeatModeler.  
Run this on AWS node configured for RM:  
    ec2-13-59-253-165.us-east-2.compute.amazonaws.com  
```
ssh  alexl@ec2-13-59-253-165.us-east-2.compute.amazonaws.com
# set the environment
umask 002
species="Fragaria_vesca_4"
cd /data
mkdir -p $species
cd $species
mkdir -p data
cd data
scp alexl@topaz.gatech.edu:/storage3/EukSpecies/$species/data/genome.fasta  .
  ## password
cd /data/$species/
mkdir run_a
cd run_a
mkdir RModeler RMasker
cp  ../../bin/run_masking_denovo.sh .
nohup ./run_masking_denovo.sh >&  loginfo &
# wait and check
cd RMasker
scp  genome.fasta.masked  alexl@topaz.gatech.edu:/storage3/EukSpecies/$species/data
  ## password
exit
```
Get masking coordinates from soft-masked sequence 
```
cd $base/annot/
soft_fasta_to_3 < ../data/genome.fasta.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3 "\t.\t.\t.\t." }' > mask.gff

# collect stat
cat mask.gff | awk '{sum+=($5-$4+1)}END{print sum}'
cat mask.gff | awk '{if ( $5-$4+1 >= 1000) sum+=($5-$4+1)}END{print sum}'

```
### Annotation  
Download annotation from 
Select only protein coding genes from annotation and save it in GFF3 and GTF (stop codon included) formats.  
```
# download
cd $base/arx
wget 

# select 
gff_to_gff_subset.pl --in Fragaria_vesca_v4.0.a1.transcripts.gff3 --out annot.gff3 --list list.tbl --col 2

# add GFF3 headers
echo "##gff-version 3" > tmp_annot.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> tmp_annot.gff3
cat annot.gff3 | grep -v gff-version  >> tmp_annot.gff3
mv tmp_annot.gff3 annot.gff3

# check
/home/tool/gt/bin/gt  gff3validator annot.gff3

# make nice
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

# separate pseudo
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3

# add features
enrich_gff.pl --in annot.gff3 --out ref.gff3 --cds --seq ../data/genome.fasta --v --warnings
gff3_to_gtf.pl ref.gff3 ref.gtf

# check
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gff3
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gtf
/home/braker/src/eval-2.2.8/validate_gtf.pl  ref.gtf

# move files to annot folder
mv ref.gff3     ../annot/annot.gff3
mv ref.gtf      ../annot/annot.gtf
mv pseudo.gff3  ../annot/

rm annot.gff3
gzip Fragaria_vesca_v4.0.a1.transcripts.gff3
```

