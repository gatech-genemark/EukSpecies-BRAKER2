# Species: _Bombus terrestris_

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

```bash
cd $base/arx
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.fna.gz
gunzip GCF_000214255.1_Bter_1.0_genomic.fna.gz

cat GCF_000214255.1_Bter_1.0_genomic.fna  | grep '^>' > defline
cat defline | cut -f1 -d ' ' | cut -b2- > z
paste z z > list.tbl

get_fasta_with_tag.pl --swap --in GCF_000214255.1_Bter_1.0_genomic.fna  --out tmp_genome.fasta  --list list.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

rm tmp_genome.fasta
mv genome.fasta ../data/genome.fasta
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

According to BeeBase, "The NCBI Annotation is the Official Gene Set of the Bumble Bee Genome Consortium." (http://hymenopteragenome.org/beebase/?q=consortium_datasets).

```bash
cd $base/arx
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.gff.gz
gunzip GCF_000214255.1_Bter_1.0_genomic.gff.gz
echo "##gff-version 3" > annot.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' |  grep -v '^$' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> annot.gff3
cat GCF_000214255.1_Bter_1.0_genomic.gff | grep -v "#" >> annot.gff3

gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3 annot.gff3

select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3
mv pseudo.gff3 ../annot

# Works when line 943 in enrich_gff.pl is commented out
enrich_gff.pl --in annot.gff3 --out tmp_annot.gff3 --cds --seq ../data/genome.fasta --v --warnings
mv tmp_annot.gff3 annot.gff3

gff3_to_gtf.pl annot.gff3 annot.gtf

# Check
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf

mv annot.gff3     ../annot
mv annot.gtf      ../annot
```

### Partial genes

```bash
cd $base/annot
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.gtf.gz
gunzip GCF_000214255.1_Bter_1.0_genomic.gtf.gz
grep "partial \"true\"" GCF_000214255.1_Bter_1.0_genomic.gtf |  grep -o -P "transcript_id\ \"[^\"]+\"" | cut -f2 -d " " | tr -d \" | sort | uniq > partial_transcripts
grep "partial \"true\"" GCF_000214255.1_Bter_1.0_genomic.gtf |  grep -o -P "gene_id\ \"[^\"]+\"" | cut -f2 -d " " | tr -d \" | sort | uniq > partial_genes

grep -Ff partial_transcripts annot.gtf > incompleteTranscripts.gtf
grep -v -Ff partial_transcripts annot.gtf > completeTranscripts.gtf
grep -v -Ff partial_genes annot.gtf > completeGenes.gtf
grep -Ff partial_genes annot.gtf > incompleteGenes.gtf
```
