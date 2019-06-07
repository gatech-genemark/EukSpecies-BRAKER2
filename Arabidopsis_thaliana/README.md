# Species: _Arabidopsis thaliana_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Environmental variables setup:  
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
```
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip  GCF_000001735*.fna.gz

grep '^>' GCF*.fna
grep '^>' GCF*.fna | sed 's/ .*//' | tr -d '>' | grep -v NC_037304 | grep -v NC_000932  > list.tbl
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
cd /data
mkdir Arabidopsis_thaliana
cd Arabidopsis_thaliana
mkdir data RModeler RMasker
cd data
scp alexl@topaz.gatech.edu:/storage3/w/alexl/EukSpecies/Arabidopsis_thaliana/data/genome.fasta .
  ## passwd
cd ../RModeler
cp ../../bin/run_RModeler.sh .
./run_RModeler.sh
# wait
# check
cd ../RMasker
ln -s ../data/genome.fasta 
run_RMasker.sh
# wait
# check
scp zzz alexl@topaz.gatech.edu:/storage3/w/alexl/EukSpecies/Arabidopsis_thaliana/data
  # passwd
```



