# Species: _Solanum lycopersicum_

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

ITAG 4.0 assembly

```bash
cd $base/arx
wget ftp://ftp.solgenomics.net/tomato_genome/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa

grep '^>' S_lycopersicum_chromosomes.4.00.fa  > deflines

cat deflines | cut -b2- > tmp_1
cat deflines | cut -b7- | sed 's/ch/chr_/' > tmp_2
paste tmp_1 tmp_2 | grep -v "00" > list.tbl
rm tmp_1 tmp_2

probuild --stat --seq  S_lycopersicum_chromosomes.4.00.fa  --details
probuild --stat_fasta --seq  S_lycopersicum_chromosomes.4.00.fa  --details

get_fasta_with_tag.pl --swap --in S_lycopersicum_chromosomes.4.00.fa  --out genome.fasta  --list list.tbl --v

probuild --stat --details --seq genome.fasta
probuild --stat_fasta --details --seq genome.fasta

probuild --reformat_fasta --in genome.fasta --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --original

probuild --stat --details --seq ../data/genome.fasta
probuild --stat_fasta --details --seq  ../data/genome.fasta

rm genome.fasta

gzip S_lycopersicum_chromosomes.4.00.fa
```

# Annotation

```bash
cd $base/arx

wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_gene_models.gff
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_REPET_repeats_aggressive.gff
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_RepeatModeler_repeats_light.gff

gff_to_gff_subset.pl  --swap  --list list.tbl  --v  --in ITAG4.0_gene_models.gff  --out itag.gff
echo "##gff-version 3" > tmp_itag.gff
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | grep chr_ | tr -d '>' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> tmp_itag.gff
cat itag.gff | grep -v '##gff-version'  >> tmp_itag.gff
mv tmp_itag.gff itag.gff

# check
gt  gff3validator itag.gff

# reformat
gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o itag_all.gff3  itag.gff
enrich_gff.pl --in itag_all.gff3 --out itag.gff3 --cds --seq ../data/genome.fasta --v --warnings
gff3_to_gtf.pl itag.gff3 itag.gtf
rm itag_all.gff3 itag.gff

mv itag.gff3 ../annot/annot.gff3
mv itag.gtf  ../annot/annot.gtf

# move masking coordinates to new seq ID and subset it
gff_to_gff_subset.pl  --swap  --list list.tbl  --v  --out $base/annot/itag_RM_softmask.gff     --in ITAG4.0_RepeatModeler_repeats_light.gff
gff_to_gff_subset.pl  --swap  --list list.tbl  --v  --out $base/annot/itag_REPET_softmask.gff  --in ITAG4.0_REPET_repeats_aggressive.gff
```

### Categorize complete and incomplete transcripts

```bash
cd $base/annot
findPartialGenes.py annot.gtf  --completeTranscripts completeTranscripts.gtf --incompleteTranscripts incompleteTranscripts.gtf --completeGenes completeGenes.gtf --incompleteGenes incompleteGenes.gtf
```

# Masking

### Community coordinates

```bash
cd $base/data

# mask genome using REPET aggressive masking coordinates
bedtools  maskfasta  -fi genome.fasta  -bed ../annot/itag_REPET_softmask.gff  -fo genome_REPET.fasta  -soft

# mask genome using RepeatModeler light masking coordinates
/home/tool/bedtools2/bin/bedtools  maskfasta  -fi genome.fasta  -bed ../annot/itag_RM_softmask.gff  -fo genome_RM.fasta  -soft
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
