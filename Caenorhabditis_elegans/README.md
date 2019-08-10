# Species: _Caenorhabditis_elegans_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Setup environment on GT cluster as:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Caenorhabditis_elegans"

export PATH="$base/bin:$PATH"
export base="$base/$species"
cd $base
if [ "$(pwd)" != "$base" ]; then echo "error, folder not found: $base"; fi
```
Create core folders  
```
cd $base
mkdir -p arx annot data mask
```
Download genomic sequence and reformat it:  
 * simplified FASTA defline with a first word in defline as a unique sequence ID
 * select only nuclear DNA (exclude organelles)
 * set sequence in all uppercase
When possible use genomic sequence from NCBI.  
Match sequence ID in FASTA file with sequence ID in annotation file.  
Use ID from annotation.  
Keep IDs in the file "list.tbl".  
First column in the table is sequence ID and second column is annotation ID.  

Description of assembly is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6  
```
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
gunzip  GCF_000002985*.fna.gz

grep '^>' GCF*.fna
grep '^>' GCF*.fna | grep -v NC_001328 | cut -f1,5 -d' ' | sed 's/^>//'  > list.tbl
get_fasta_with_tag.pl --swap --in GCF_000002985.6_WBcel235_genomic.fna  --out genome.fasta  --list list.tbl --v
probuild --stat --details --seq genome.fasta
probuild --reformat_fasta --in genome.fasta --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --original
rm genome.fasta
probuild --stat --details --seq ../data/genome.fasta

gzip  GCF_000002985*.fna
```
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
cd ..
cp ../bin/run_masking.sh .
nohup ./run_masking.sh >&  loginfo &
# wait and check
cd RMasker
scp  genome.fasta.masked  alexl@topaz.gatech.edu:/storage3/w/alexl/EukSpecies/$species/data
  ## password
exit
```
Download annotation from WormBase.  
NCBI RefSeq is using annotation from WormBase.  
Select only protein coding genes from annotation and save it in GFF3 and GTF (stop codon included) formats.  
```
cd $base/arx
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/gff/c_elegans.PRJNA13758.WS271.annotations.gff3.gz
gunzip  c_elegans.PRJNA13758.WS271.annotations.gff3.gz

gff_to_gff_subset.pl  --in c_elegans.PRJNA13758.WS271.annotations.gff3  --out tmp_annot.gff3 --list list.tbl --col 2
grep '^#' tmp_annot.gff3  | sort -k2,2 | uniq > annot.gff3
cat tmp_annot.gff3 | grep -P '\tWormBase\t' >> annot.gff3
rm tmp_annot.gff3
#check
/home/tool/gt/bin/gt  gff3validator annot.gff3
#reformat
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

enrich_gff.pl --in annot.gff3 --out ref.gff3 --cds 
gff3_to_gtf.pl ref.gff3 ref.gtf
ln -s ref.gtf annot.gtf

# check
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf
/home/braker/src/eval-2.2.8/validate_gtf.pl -c annot.gtf

mv annot.gff3 ../annot/
mv annot.gtf  ../annot/
mv ref.gff3   ../annot/
mv ref.gtf    ../annot

gzip c_elegans.PRJNA13758.WS271.annotations.gff3

# separate pseudo
cd $base/annot/
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3

# masking coordinates
cd $base/annot/
soft_fasta_to_3 < ../data/genome.fasta.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3+1 "\t.\t.\t.\t." }' > mask.gff

```


