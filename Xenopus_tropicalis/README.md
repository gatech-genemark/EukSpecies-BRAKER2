# Species: _Xenopus tropicalis_

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
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz
gunzip GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz

grep '^>' GCF_000004195.4_UCB_Xtro_10.0_genomic.fna  > deflines
cat deflines | grep -v "^>NW" | grep -v "^>NC_006839.1" | cut -b2- | awk '{print $1 "\t" $7}' | tr -d ',' > list.tbl

../../bin/get_fasta_with_tag.pl --swap --in GCF_000004195.4_UCB_Xtro_10.0_genomic.fna  --out tmp_genome.fasta  --list list.tbl --v
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

Additional masking by Tandem Repeats Finder (with maximum repeat period size = 500) was applied
since the default run of RepeatMasker/RepeatModeler did not identify a significant portion of long
tandem repeats (with repeat pattern length > 10) in this genome.

[Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) (TRF) v4.07b was run with the following command:

```bash
cd $base/data
trf genome.fasta 2 7 7 80 10 50 500 -d -m -h
```

**TODO** add link to the script

Coordinates of repeats from TRF .dat format were converted to .gff with a custom [parseTrfOutput.py]() script

```bash
parseTrfOutput.py genome.fasta.2.7.7.80.10.50.500.dat --minCopies 1 --statistics STATS \
    > genome.fasta.2.7.7.80.10.50.500.raw.gff
```

The masking coordinates were applied on top of RepeatModeler/RepeatMasker masking

```bash
# Sort gff
sort -k1,1 -k4,4n -k5,5n genome.fasta.2.7.7.80.10.50.500.raw.gff > sorted
# Merge overlapping repeats
bedtools merge -i sorted | awk ’BEGIN{OFS="\t"} {print $1,"trf","repeat",$2+1,$3,".",".",".","."}’ \
    > genome.fasta.2.7.7.80.10.50.500.merged.gff
# Apply masking
bedtools maskfasta -fi genome.fasta.masked -bed genome.fasta.2.7.7.80.10.50.500.merged.gff \
    -fo genome.fasta.combined.masked -soft

cd $base/annot
soft_fasta_to_3 < ../data/genome.fasta.combined.masked | awk '{print $1 "\tsoft_masking\trepeat\t" $2+1 "\t" $3 "\t.\t.\t.\t." }' > combined.mask.gff
```

# Annotation

NCBI

```bash
cd $base/arx
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_genomic.gff.gz
gunzip GCF_000004195.4_UCB_Xtro_10.0_genomic.gff.gz

gff_to_gff_subset.pl  --in GCF_000004195.4_UCB_Xtro_10.0_genomic.gff  --out tmp_annot.gff3  --list list.tbl  --col 1  --v --swap
echo "##gff-version 3" > annot.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' |  grep -v '^$' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> annot.gff3
cat tmp_annot.gff3 | grep -v "##" >> annot.gff3

# Remove miRNA genes because they do not have valid parents in gff3
grep -v gene=mir.* annot.gff3 > tmp_annot.gff3; mv tmp_annot.gff3 annot.gff3

# Remove some other invalid entries
grep -v exon-XR_004221006.1-1 annot.gff3 > tmp_annot.gff3; mv tmp_annot.gff3 annot.gff3
grep -v rna-XR_004221000.1 annot.gff3 > tmp_annot.gff3; mv tmp_annot.gff3 annot.gff3
grep -v "small Cajal body-specific RNA" annot.gff3 > tmp_annot.gff3; mv tmp_annot.gff3 annot.gff3

# Validate
gt gff3validator annot.gff3

gt gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_annot.gff3  annot.gff3
mv tmp_annot.gff3 annot.gff3

select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3
mv pseudo.gff3 ../annot

# Works when line 943 in enrich_gff.pl is commented out
enrich_gff.pl --in annot.gff3 --out tmp_annot.gff3 --cds
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
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_genomic.gtf.gz
gunzip GCF_000004195.4_UCB_Xtro_10.0_genomic.gtf.gz

grep "partial \"true\"" GCF_000004195.4_UCB_Xtro_10.0_genomic.gtf |  grep -o -P "transcript_id\ \"[^\"]+\"" | cut -f2 -d " " | tr -d \" | sort | uniq > partial_transcripts
grep "partial \"true\"" GCF_000004195.4_UCB_Xtro_10.0_genomic.gtf |  grep -o -P "gene_id\ \"[^\"]+\"" | cut -f2 -d " "  | tr -d \" |awk '{print "-" $1 "\""}' > partial_genes

grep -Ff partial_transcripts annot.gtf > incompleteTranscripts.gtf
grep -v -Ff partial_transcripts annot.gtf > completeTranscripts.gtf
grep -v -Ff partial_genes annot.gtf > completeGenes.gtf
grep -Ff partial_genes annot.gtf > incompleteGenes.gtf

cp annot.gtf annot_raw.gtf
```
