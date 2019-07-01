# Species: _Drosophila_melanogaster_  
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup  
Project is set in bash shell.  

Setup environment on GT cluster as:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Drosophila_melanogaster"

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

Description of assembly is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4  
```
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
gunzip  GCF_000001215*.fna.gz

grep '^>' GCF*.fna
grep '^>' GCF*.fna | grep -v NW_00  | grep -v NC_024511 | cut -f1,5 -d' ' | sed 's/^>//' > list.tbl
get_fasta_with_tag.pl  --swap  --in GCF_000001215*.fna  --out tmp_genome.fasta  --list list.tbl  --v
probuild --stat --details --seq tmp_genome.fasta
probuild --reformat_fasta --in tmp_genome.fasta --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --original
rm tmp_genome.fasta
probuild --stat --details --seq ../data/genome.fasta

gzip  GCF_000001215*.fna
```
Run _de novo_ masking of genome using RepeatModeler.  
Run this on AWS node configured for RM:  
    ec2-13-59-253-165.us-east-2.compute.amazonaws.com
```
ssh  alexl@ec2-13-59-253-165.us-east-2.compute.amazonaws.com
# set the environment
umask 002
species="Drosophila_melanogaster"

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
Download annotation from FlyBase.  
NCBI RefSeq is using annotation from FlyBase.  
Select protein coding genes from annotation and save them in GFF3 and GTF (stop codon included) formats.  
```
cd $base/arx
wget ftp://ftp.flybase.net/releases/FB2019_03/dmel_r6.28/gff/dmel-all-no-analysis-r6.28.gff.gz
gunzip  dmel-all-no-analysis-*.gff.gz

gff_to_gff_subset.pl  --in dmel-all-no-analysis-r6.28.gff  --out tmp.gff3  --list list.tbl  --col 2
#check
/home/tool/gt/bin/gt  gff3validator  tmp.gff3
# some CDS lines with multiple parents require different phases : this creates a new CDS line with corrected phase
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o annot.gff3  tmp.gff3
compare_intervals_exact.pl  --f1 annot.gff3  --f2 tmp.gff3
rm  tmp.gff3

gff3_to_gtf.pl  annot.gff3  annot.gtf
# eval validate is not handling correctly transspliced genes on different strands
# do not use "-f" fixed by eval version
/home/braker/src/eval-2.2.8/validate_gtf.pl -c annot.gtf
compare_intervals_exact.pl  --f1 annot.gff3  --f2 annot.gtf

mv annot.gff3 ../annot/
mv annot.gtf  ../annot/

gzip dmel-all-no-analysis-*.gff
```


