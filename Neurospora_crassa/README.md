# Species: _Neurospora crassa_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup  
```
species="Neurospora_crassa"
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
Assembly description is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000182925.2  
```
# download data
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz
gunzip  GCF_000182925*.fna.gz

# create ID table
grep '^>' GCF*.fna > deflines
cat deflines | grep -v '>NW_' | grep -v NC_026614 | cut -f1,7 -d' ' | tr -d ',' | tr -d '>' | sed 's/ / Supercontig_ /' > list.tbl

# Modify manually sequence ID's in the second column of list.tbl to match the annoatation file

# select and reformat sequence; all uppercase 
get_fasta_with_tag.pl --swap --in GCF_000182925.2_NC12_genomic.fna  --out tmp_genome.fasta  --list list.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

# check sequence and clean folder
probuild --stat --details --seq tmp_genome.fasta
probuild --stat --details --seq genome.fasta

# put in work folder
mv genome.fasta ../data/genome.fasta

# clean tmp files
rm tmp_genome.fasta
gzip  GCF_000182925*.fna
```
### Masking: _de novo_ and _species specific_
Run _de novo_ masking of genome using RepeatModeler.  
Run this on AWS node configured for RM:  
    ec2-13-59-253-165.us-east-2.compute.amazonaws.com  
```
ssh  alexl@ec2-13-59-253-165.us-east-2.compute.amazonaws.com
# set the environment
umask 002
species="Neurospora_crassa"

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
Get masking coordinates from soft-masked sequence 
```
cd $base/annot/
soft_fasta_to_3 < ../data/genome.fasta.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3+1 "\t.\t.\t.\t." }' > mask.gff
```
### Annotation  
Download annotation from JGI.  
JGI annotation is copy of Broad annotation.
https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Neucr2  
NCBI RefSeq may differ from Broad annotation.  
Select only protein coding genes from annotation and save in GFF3 and GTF (stop codon included) formats.  
```
# manual download from JGI website above
cd $base/arx
gunzip  Neucr2.filtered_proteins.BroadModels.gff3.gz

# select
gff_to_gff_subset.pl  --in Neucr2.filtered_proteins.BroadModels.gff3  --out annot.gff3 --list list.tbl --col 2

# reformat into "nice" gff3
grep '^#' annot.gff3 > tmp_annot.gff3
grep -v '^#' annot.gff3 >> tmp_annot.gff3
mv tmp_annot.gff3 annot.gff3

# check
/home/tool/gt/bin/gt  gff3validator annot.gff3

# make nice
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

# separate pseudo
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3

enrich_gff.pl --in annot.gff3 --out ref.gff3 --cds --seq ../data/genome.fasta --v --warnings
gff3_to_gtf.pl ref.gff3 ref.gtf

# check
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gff3
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gtf
/home/braker/src/eval-2.2.8/validate_gtf.pl  ref.gtf

mv ref.gff3   ../annot/annot.gff3
mv ref.gtf    ../annot/annot.gtf
mv pseudo.gff3  ../annot/

rm annot.gff3
gzip Neucr2.filtered_proteins.BroadModels.gff3
```
### Create a version of annotation with NCBI contig names
```
cd $base/annot/
cp annot.gtf annot_NC_contigs.gtf
for i in {1..7}; do sed -i "s/Supercontig_$i/NC_02650${i}.1/" annot_NC_contigs.gtf; done
```
