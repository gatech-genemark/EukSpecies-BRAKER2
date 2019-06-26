# Species: _Caenorhabditis_elegans_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Environmental variables setup on GT cluster:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Caenorhabditis_elegans"

export PATH="$base/bin:$PATH"
export base="$base/$species"
cd $base
if [ "$(pwd)" != "$base" ]; then echo "error, folder not found: $base"; fi
```
Create core folder structure
```
cd $base
mkdir -p arx annot data
```
Download genomic sequence and reformat it:  
 * simplified FASTA defline with unique sequence ID as a first word in defline
 * select only nuclear DNA (exclude organelles)
 * all uppercase

Use genomic sequence from NCBI, when possible.  
Match sequence ID's in FASTA file with sequence ID's in annotation file.  
Use ID's sequence from annotation.  
Assembly description is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6
Keep FASTA IDs in the file "list.tbl".  
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

gff_to_gff_subset.pl  --in c_elegans.PRJNA13758.WS271.annotations.gff3  --out tmp.gff3 --list list.tbl --col 2
grep '^#' tmp.gff3  | sort -k2,2 | uniq > annot.gff3
cat tmp.gff3 | grep -P '\tWormBase\t' >> annot.gff3
rm tmp.gff3
#check
/home/tool/gt/bin/gt  gff3validator annot.gff3
#reformat
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  annot.gff3
mv tmp.gff3  annot.gff3

gff3_to_gtf.pl annot.gff3 annot.gtf
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf
#check
/home/braker/src/eval-2.2.8/validate_gtf.pl -c annot.gtf

mv annot.gff3 ../annot/
mv annot.gtf  ../annot/

gzip c_elegans.PRJNA13758.WS271.annotations.gff3
```


