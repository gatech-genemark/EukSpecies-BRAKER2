# Species: _Populus trichocarpa_

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
# Download instructions, Alex, please fill in

grep ">" Ptrichocarpa_533_v4.0.softmasked.fa > deflines
cat deflines | grep Chr | tr -d ">" | awk '{print $1 "\t" $1}' > list.tbl

get_fasta_with_tag.pl --swap --in Ptrichocarpa_533_v4.0.softmasked.fa  --out tmp_genome.fasta  --list list.tbl --v
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

TODO
