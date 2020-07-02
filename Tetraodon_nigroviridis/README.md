# Species: _Tetraodon nigroviridis_

Tomas Bruna, Alex Lomsadze
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
wget ftp://ftp.ensembl.org/pub/release-99/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.dna.toplevel.fa.gz
gunzip Tetraodon_nigroviridis.TETRAODON8.dna.toplevel.fa.gz

grep '^>' Tetraodon_nigroviridis.TETRAODON8.dna.toplevel.fa  > deflines
cat deflines  | cut -f1 -d' ' | grep -v MT | cut -b2- > z
paste z z > list.tbl
rm z

get_fasta_with_tag.pl --swap --in Tetraodon_nigroviridis.TETRAODON8.dna.toplevel.fa   --out tmp_genome.fasta  --list list.tbl --v
probuild --reformat_fasta --in tmp_genome.fasta --out genome.fasta --uppercase 1 --letters_per_line 60 --original

rm tmp_genome.fasta
mv genome.fasta  ../data
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
cd $base/arx
wget ftp://ftp.ensembl.org/pub/release-99/gff3/tetraodon_nigroviridis/Tetraodon_nigroviridis.TETRAODON8.99.gff3.gz
gunzip Tetraodon_nigroviridis.TETRAODON8.99.gff3.gz

gff_to_gff_subset.pl  --in Tetraodon_nigroviridis.TETRAODON8.99.gff3  --out tmp_annot.gff3  --list list.tbl --col 2 --v
echo "##gff-version 3" > annot.gff3
probuild --stat_fasta --seq ../data/genome.fasta | cut -f1,2 | tr -d '>' | grep -v '^$' | awk '{print "##sequence-region  " $1 "  1 " $2}' >> annot.gff3
cat tmp_annot.gff3 | grep -v "^#" >> annot.gff3
rm  tmp_annot.gff3

# check gff3 validity
gt  gff3validator annot.gff3

# reformat
gt gff3 -force -tidy -sort -retainids -checkids -o tmp_annot.gff3 annot.gff3
mv tmp_annot.gff3  annot.gff3

select_pseudo_from_nice_gff3.pl annot.gff3 pseudo.gff3
mv pseudo.gff3 ../annot

# enrich
enrich_gff.pl --in annot.gff3 --out tmp_annot.gff3 --cds --seq ../data/genome.fasta --v --warnings
mv tmp_annot.gff3  annot.gff3

# to gtf
gff3_to_gtf.pl annot.gff3 annot.gtf

# checks
compare_intervals_exact.pl --f1 annot.gff3  --f2 annot.gtf

mv annot.gff3     ../annot
mv annot.gtf      ../annot
```

### Categorize complete and incomplete transcripts

```bash
cd $base/annot
findPartialGenes.py annot.gtf  --completeTranscripts completeTranscripts.gtf --incompleteTranscripts incompleteTranscripts.gtf --completeGenes completeGenes.gtf --incompleteGenes incompleteGenes.gtf
```

### Dealing with incomplete CDS

**This is an old description of an alternative way to find partial genes, left here for historical reasons.**

* The following script:
    * Flags partial CDS
    * Removes extra start and stops in the enriched annotation
    * Splits the annotation into files with:
        * Complete/incomplete transcripts
        * Complete/incomplete genes. Gene is considered to be incomplete if at least one of its transcripts is incomplete.

Assumes that annot.gtf is the enriched version of annotation with **incorrect starts and stops being part of partial CDS segments**.

```bash
cd $base/arx
wget ftp://ftp.ensembl.org/pub/release-99/gtf/tetraodon_nigroviridis/Tetraodon_nigroviridis.TETRAODON8.99.gtf.gz
gunzip Tetraodon_nigroviridis.TETRAODON8.99.gtf.gz

cd $base/annot
flagPartialCDSFromEnsembl.py ../arx/Tetraodon_nigroviridis.TETRAODON8.99.gtf annot.gtf --incompleteTranscriptsOutput incompleteTranscripts.gtf \
    --completeTranscriptsOutput completeTranscripts.gtf --fullOutput annot_fixed_partial.gtf --completeGenesOutput completeGenes.gtf \
    --incompleteGenesOutput incompleteGenes.gtf
mv annot.gtf annot_raw.gtf
mv annot_fixed_partial.gtf annot.gtf
```
