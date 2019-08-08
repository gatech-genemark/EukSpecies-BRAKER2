# Species: _Arabidopsis thaliana_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Setup environment on GT cluster as:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Arabidopsis_thaliana"

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

Description of assembly is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.4  
```
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip  GCF_000001735*.fna.gz

grep '^>' GCF*.fna
grep '^>' GCF*.fna  | grep -v NC_000932 | grep -v NC_037304 | cut -f1,5 -d' ' | sed 's/^>//'  | sed 's/ / Chr/' > list.tbl
get_fasta_with_tag.pl --swap --in GCF_000001735.4_TAIR10.1_genomic.fna  --out genome.fasta  --list list.tbl --v
probuild --stat --details --seq genome.fasta
probuild --reformat_fasta --in genome.fasta --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --original
rm genome.fasta
probuild --stat --details --seq ../data/genome.fasta

gzip  GCF_000001735*.fna
```
Run _de novo_ masking of genome using RepeatModeler.  
Run this on AWS node configured for RM:  
    ec2-13-59-253-165.us-east-2.compute.amazonaws.com
```
ssh  alexl@ec2-13-59-253-165.us-east-2.compute.amazonaws.com
# set the environment
umask 002
species="Arabidopsis_thaliana"

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
Download annotation from Tair.  
NCBI RefSeq is using annotation from Tair.  
Select only protein coding genes from annotation and save it in GFF3 and GTF (stop codon included) formats.  
```
cd $base/arx
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
gunzip  Araport11_GFF3_genes_transposons.201606.gff.gz

gff_to_gff_subset.pl  --in Araport11_GFF3_genes_transposons.201606.gff  --out tmp_annot.gff3 --list list.tbl --col 2
echo "##gff-version 3" > annot.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' | grep Chr | awk '{print "##sequence-region  " $1 "  1 " $2}' >> annot.gff3
cat tmp_annot.gff3 | grep -v gff-version  >> annot.gff3
rm  tmp_annot.gff3
#check
/home/tool/gt/bin/gt  gff3validator annot.gff3
# reformat
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

gff3_to_gtf.pl annot.gff3 annot.gtf
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf
# check
/home/braker/src/eval-2.2.8/validate_gtf.pl -c annot.gtf

mv annot.gff3 ../annot/
mv annot.gtf  ../annot/

gzip Araport11_GFF3_genes_transposons.201606.gff

# separate pseudo
cd $base/annot/
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3

# masking coordinates
cd $base/annot/
soft_fasta_to_3 < ../data/genome.fasta.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3+1 "\t.\t.\t.\t." }' > mask.gff
```


