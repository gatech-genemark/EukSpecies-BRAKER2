# Species: _Caenorhabditis_elegans_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup  
```
species="Caenorhabditis_elegans"
base="/storage3/w/alexl/EukSpecies"
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
Assembly description is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6  
```
# download data
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
gunzip  GCF_000002985*.fna.gz

# create ID table
grep '^>' GCF*.fna > deflines
cat deflines | grep -v NC_001328 | cut -f1,5 -d' ' | sed 's/^>//'  > list.tbl

# select and reformat sequence
get_fasta_with_tag.pl --swap --in GCF_000002985.6_WBcel235_genomic.fna  --out tmp_genome.fasta  --list list.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

# check sequence and clean folder
probuild --stat --details --seq tmp_genome.fasta
probuild --stat --details --seq genome.fasta

# put in work folder
mv genome.fasta ../data/genome.fasta

# clean tmp files
rm defline tmp_genome.fasta
gzip  GCF_000001735*.fna
```
### Masking: _de novo_ and _species specific_
Run _de novo_ masking of genome using RepeatModeler.  
Run this on AWS node configured for RM:  
    ec2-13-59-253-165.us-east-2.compute.amazonaws.com  
```
ssh  alexl@ec2-13-59-253-165.us-east-2.compute.amazonaws.com
# set the environment
umask 002
species="Caenorhabditis_elegans"

cd /data
mkdir -p $species
cd $species
mkdir -p data RModeler RMasker
cd data
scp alexl@topaz.gatech.edu:/storage3/w/alexl/EukSpecies/$species/data/genome.fasta  .
  ## password
cd /data/$species/
cp ../bin/run_masking.sh .
nohup ./run_masking.sh >&  loginfo &
# wait and check
cd RMasker
scp  genome.fasta.masked  alexl@topaz.gatech.edu:/storage3/w/alexl/EukSpecies/$species/data
  ## password
exit
```
### Annotation 
Download annotation from WormBase. NCBI RefSeq is using annotation from WormBase.  
Select only protein coding genes from annotation and save it in GFF3 and GTF (stop codon included) formats.  
```
# download
cd $base/arx
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS271.annotations.gff3.gz
gunzip  c_elegans.PRJNA13758.WS271.annotations.gff3.gz

# select 
gff_to_gff_subset.pl  --in c_elegans.PRJNA13758.WS271.annotations.gff3  --out annot.gff3 --list list.tbl --col 2

# reformat into "nice" gff3
echo "##gff-version 3" > tmp_annot.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' | grep -v '^$' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> tmp_annot.gff3
cat annot.gff3 | grep -P '\tWormBase\t' >> tmp_annot.gff3
mv tmp_annot.gff3 annot.gff3

# check
/home/tool/gt/bin/gt  gff3validator annot.gff3

# make nice
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

enrich_gff.pl --in annot.gff3 --out ref.gff3 --cds 
gff3_to_gtf.pl ref.gff3 ref.gtf

# check
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gff3
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gtf
/home/braker/src/eval-2.2.8/validate_gtf.pl  ref.gtf

mv ref.gff3   ../annot/annot.gff3
mv ref.gtf    ../annot/annot.gtf

gzip c_elegans.PRJNA13758.WS271.annotations.gff3

# separate pseudo
cd $base/annot/
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3

# masking coordinates
cd $base/annot/
soft_fasta_to_3 < ../data/genome.fasta.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3+1 "\t.\t.\t.\t." }' > mask.gff
```


