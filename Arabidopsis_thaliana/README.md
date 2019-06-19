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
```
Run _de novo_ masking of genome using RepeatMOdeler.  
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



