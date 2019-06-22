# Species: _Arabidopsis thaliana_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Environmental variables setup on GT cluster:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Arabidopsis_thaliana"

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
Assembly description is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.4
Keep FASTA IDs in the file "list.tbl".  
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
cd ../RModeler
cp ../../bin/run_RModeler.sh .
./run_RModeler.sh
# wait and check
cd ../RMasker
ln -s ../data/genome.fasta
cp ../../bin/run_RMasker.sh .
./run_RMasker.sh
# wait and check
scp  genome.fasta.masked  alexl@topaz.gatech.edu:/storage3/w/alexl/EukSpecies/$species/data
  ## password
```
Download annotation from Tair.  
NCBI RefSeq is using annotation from Tair.  
Select only protein coding genes from annotation and save it in GFF3 and GTF (stop codon included) formats.  
```
cd $base/arx
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
gunzip  Araport11_GFF3_genes_transposons.201606.gff.gz

gff_to_gff_subset.pl  --in Araport11_GFF3_genes_transposons.201606.gff  --out annot.gff3 --list list.tbl --col 2
echo "##gff-version 3" > tmp_header.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' | grep Chr | awk '{print "##sequence-region  " $1 "  1 " $2}' >> tmp_header.gff3
cat annot.gff3 | grep -v gff-version  >> tmp_header.gff3
mv tmp_header.gff3  annot.gff3
gff3_to_gtf.pl annot.gff3 annot.gtf
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf

# optional testing and reformating
/home/tool/gt/bin/gt  gff3validator annot.gff3
/home/tool/gt/bin/gt  gff3  -addintrons  -sort  -checkids  -o test.gff3  -retainids  -force  annot.gff3
mv test.gff3 annot.gff3
/home/braker/src/eval-2.2.8/validate_gtf.pl -c -f annot.gtf
mv annot.fixed.gtf annot.gtf
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf

mv annot.gff3 ../annot/
mv annot.gtf  ../annot/

gzip Araport11_GFF3_genes_transposons.201606.gff
```


