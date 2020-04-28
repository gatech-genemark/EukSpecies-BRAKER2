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

TODO

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
