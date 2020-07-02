# Species: _Danio rerio_

Alex Lomsadze, Tomas Bruna,
Georgia Institute of Technology,
2020

## Project setup

```bash
base=$(pwd)
export PATH="$base/../bin:$PATH"
umask 002
# Create core folders
mkdir arx annot data
```

# Genome

Description of assembly is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.6/  

```bash
cd $base/arx
mkdir ensembl refseq

cd $base/arx/refseq
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz
gunzip GCF_000002035.6_GRCz11_genomic.fna.gz

grep '^>' GCF*.fna > deflines.refseq
cat deflines.refseq | grep -v "^>NW" | grep -v "^>NC_002333.2" | cut -b2- | awk '{print $1 "\t" $7}' | tr -d ',' > list_refseq.tbl

get_fasta_with_tag.pl --swap --in GCF_000002035.6_GRCz11_genomic.fna  --out tmp_genome.fasta  --list list_refseq.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

rm tmp_genome.fasta
mv genome.fasta ../../data/genome.fasta
gzip GCF_000002035.6_GRCz11_genomic.fna
```

### _De novo_ Masking

Mask genome with RepeatModeler (v open-1.0.11) and RepeatMasker (v 1.332)

```bash
cd $base/data
BuildDatabase -engine wublast -name genome genome.fasta
RepeatModeler -engine wublast -database genome
RepeatMasker -engine wublast -lib genome-families.fa -xsmall genome.fasta
```

Get masking coordinates from soft-masked sequence

```bash
cd $base/annot
soft_fasta_to_3 < ../data/genome.fasta.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3 "\t.\t.\t.\t." }' > mask.gff

# collect stat
cat mask.gff | awk '{sum+=($5-$4+1)}END{print sum}'
cat mask.gff | awk '{if ( $5-$4+1 >= 1000) sum+=($5-$4+1)}END{print sum}'
```

The masking coordinates can be applied to the unmasked genome with the following command:

```bash
cd $base/data
bedtools maskfasta -fi genome.fasta -bed ../annot/mask.gff -fo genome.fasta.masked -soft
```

# Annotation

```bash
cd $base/arx/ensembl
wget ftp://ftp.ensembl.org/pub/release-97/gff3/danio_rerio/Danio_rerio.GRCz11.97.gff3.gz
gunzip Danio_rerio.GRCz11.97.gff3.gz

wget ftp://ftp.ensembl.org/pub/release-97/gtf/danio_rerio/Danio_rerio.GRCz11.97.gtf.gz
gunzip Danio_rerio.GRCz11.97.gtf.gz

gff_to_gff_subset.pl  --in Danio_rerio.GRCz11.97.gff3  --out tmp_annot.gff3  --list ../refseq/list_refseq.tbl  --col 2  --v
echo "##gff-version 3" > annot.gff3
probuild --stat_fasta --seq ../../data/genome.fasta | cut -f1,2 | tr -d '>' | grep -v '^$' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> annot.gff3
cat tmp_annot.gff3 | grep -v "^#"  >> annot.gff3
rm  tmp_annot.gff3

# check
gt  gff3validator annot.gff3
# reformat
gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3  annot.gff3

# separate pseudo
select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3
mv pseudo.gff3 ../../annot/

enrich_gff.pl --in annot.gff3 --out ensembl.gff3 --cds --seq ../../data/genome.fasta --v --warnings

gff3_to_gtf.pl ensembl.gff3  ensembl.gtf

# check
compare_intervals_exact.pl --f1 ensembl.gff3  --f2 ensembl.gtf

mv ensembl.gff3   ../../annot/annot.gff3
mv ensembl.gtf    ../../annot/annot.gtf
```

###  APPRIS

Data from http://appris.bioinfo.cnio.es

```bash
cd $base/arx
# download
wget http://apprisws.bioinfo.cnio.es/pub/releases/2019_07.v29/datafiles/danio_rerio/GRCz10/appris_data.appris.txt

# get PRINCIPAL transcript ID's
cat ../annot/annot.gtf | grep -E -o 'gene:\w+' | sort | uniq  | sed s'/^gene://' > tmp_gene_names
fgrep -f tmp_gene_names  appris_data.appris.txt | grep PRINCIPAL | cut -f3 | sed 's/^/transcript:/' > appris.tbl
rm tmp_gene_names
select_by_trascript_id_from_gtf.pl  appris.tbl  ../annot/annot.gtf  appris.gtf
rm appris.tbl

# If multiple PRINCICAL transcripts are annotated per gene, then select the longest 
# In case of equal length, select the first one
get_longest_cds_gene_set.pl --in appris.gtf  --out appris.tbl -v
select_by_trascript_id_from_gtf.pl  appris.tbl  ../annot/annot.gtf  appris.gtf
rm appris.tbl

mv appris.gtf ../annot/
```

### Categorize complete and incomplete transcripts

```bash
cd $base/annot
findPartialGenes.py annot.gtf  --completeTranscripts completeTranscripts.gtf --incompleteTranscripts incompleteTranscripts.gtf --completeGenes completeGenes.gtf --incompleteGenes incompleteGenes.gtf
```

### Incomplete CDS

**This is an old description of an alternative way to find partial genes, left here for historical reasons.**

* The following script:
    * Flags partial CDS
    * Removes extra start and stops in the enriched annotation
    * Splits the annotation into files with:
        * Complete/incomplete transcripts
        * Complete/incomplete genes. Gene is considered to be incomplete if at least one of its transcripts is incomplete.

Assumes that annot.gtf is the enriched version of annotation with **incorrect starts and stops being part of partial CDS segments**.

```bash
cd $base/annot
flagPartialCDSFromEnsembl.py ../arx/ensembl/Danio_rerio.GRCz11.97.gtf annot.gtf --incompleteTranscriptsOutput incompleteTranscripts.gtf \
    --completeTranscriptsOutput completeTranscripts.gtf --fullOutput annot_fixed_partial.gtf --completeGenesOutput completeGenes.gtf \
    --incompleteGenesOutput incompleteGenes.gtf
mv annot.gtf annot_raw.gtf
mv annot_fixed_partial.gtf annot.gtf
```

Select complete genes in APPRIS.

```bash
flagPartialCDS.py ../arx/ensembl/Danio_rerio.GRCz11.97.gtf appris.gtf --fullOutput appris_fixed_partial.gtf --completeGenesOutput \
    appris_completeGenes.gtf --incompleteTranscriptsOutput /dev/null --completeTranscriptsOutput /dev/null --incompleteGenesOutput /dev/null
```
