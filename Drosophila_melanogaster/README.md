# Species: _Drosophila_melanogaster_  
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup  
```
species="Drosophila_melanogaster"
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
Assembly description is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4  
GenBank, RefSeq and FlyBase nuclear DNA sequences are identical.  
```
# download data
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
gunzip  GCF_000001215*.fna.gz

# create ID table
grep '^>' GCF*.fna > deflines
# move sequence IDs to "2L" style
cat deflines | grep -Ev 'NW_00|NC_024511'| cut -f1,5 -d' ' | sed 's/^>//' > list.tbl

# select and reformat sequence; all uppercase 
get_fasta_with_tag.pl --swap --in GCF_000001215*.fna   --out tmp_genome.fasta  --list list.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

# check sequence and clean folder
probuild --stat --details --seq tmp_genome.fasta
probuild --stat --details --seq genome.fasta

# put in work folder
mv genome.fasta ../data/genome.fasta

# clean tmp files
rm tmp_genome.fasta
gzip  GCF_000001215*.fna
```
### Masking: _de novo_ and _species specific_
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
Download annotation from FlyBase. NCBI RefSeq is using annotation from FlyBase.  
Select protein coding genes from annotation and save them in GFF3 and GTF (stop codon included) formats.  
```
# download
cd $base/arx
wget ftp://ftp.flybase.net/releases/FB2019_03/dmel_r6.28/gff/dmel-all-no-analysis-r6.28.gff.gz
gunzip  dmel-all-no-analysis-*.gff.gz

# select 
gff_to_gff_subset.pl  --in dmel-all-no-analysis-r6.28.gff  --out annot.gff3  --list list.tbl  --col 2

# check: fails on wrong format
/home/tool/gt/bin/gt  gff3validator annot.gff3

# make nice
# some CDS lines with multiple parents require different phases and thus have "." as a phase
# this splits such parents and adds phase
# some 33 CDS phases were also modifyed by this program
/home/tool/gt/bin/gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

# separate pseudo
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3

enrich_gff.pl --in annot.gff3 --out ref.gff3 --cds --seq ../data/genome.fasta --v --warnings
gff3_to_gtf.pl ref.gff3 ref.gtf

# check
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gff3
compare_intervals_exact.pl --f1 annot.gff3  --f2 ref.gtf
# eval validate is not handling correctly transspliced genes on different strands
# do not use "-f" fixed by eval version
/home/braker/src/eval-2.2.8/validate_gtf.pl  ref.gtf

# move files to annot folder
mv ref.gff3     ../annot/annot.gff3
mv ref.gtf      ../annot/annot.gtf
mv pseudo.gff3  ../annot/

rm annot.gff3
gzip dmel-all-no-analysis-*.gff
```
###  APPRIS
Data from http://appris.bioinfo.cnio.es
```
# download
wget http://apprisws.bioinfo.cnio.es/pub/releases/2019_07.v29/datafiles/drosophila_melanogaster/BDGP6/appris_data.appris.txt

# using the latest annotation - not the one used in APPRIS
# this creates small mismatch between APPRIS original and one used in this project

# get gene ID's from main annotation
cat ../annot/annot.gtf | grep -E -o 'FBgn[0-9]+' | sort | uniq > tmp_gene_names

# for each above gene ID get PRINSIPAL transcript ID
fgrep -f tmp_gene_names  appris_data.appris.txt | grep PRINCIPAL | cut -f3  > appris.tbl

# get annotation of PRINCIPAL isoforms from full annotation
# some ID's 
select_by_trascript_id_from_gtf.pl  appris.tbl  ../annot/annot.gtf  appris.gtf

rm tmp_gene_names
rm appris.tbl

# when multiple PRINCICAL transcripts are annotated per gene 
#  * select the longest
#  * in case of equal length, select the first one
get_longest_cds_gene_set.pl --in appris.gtf  --out appris_2.tbl -v
select_by_trascript_id_from_gtf.pl  appris_2.tbl  ../annot/annot.gtf  appris.gtf
rm appris_2.tbl

mv appris.gtf ../annot/

gzip appris_data.appris.txt
```
### annotation without Y chromosome
```
cd $base/annot/
grep -v "^Y" annot.gtf > annot_no_Y.gtf
grep -v "^Y" appris.gtf > appris_no_Y.gtf
```
