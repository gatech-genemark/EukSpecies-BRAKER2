# Species: _Drosophila_melanogaster_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Environmental variables setup on GT cluster: 
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Drosophila_melanogaster"

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
wget  ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
gunzip  GCF_000001215*.fna.gz

grep '^>' GCF*.fna
# adjust ID filtering
grep '^>' GCF*.fna | grep -v NW_00  | grep -v NC_024511 |  sed 's/ .*//' | tr -d '>' > list.tbl
get_fasta_with_tag.pl  --swap  --in GCF_000001215*.fna  --out genome.fasta  --list list.tbl  --v
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
species="Drosophila_melanogaster"

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

