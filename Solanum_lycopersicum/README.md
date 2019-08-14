# Species: _Solanum lycopersicum_ Tomato
Moving tomato.git to EukSpecies.git  

Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Setup environment on GT cluster as:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Solanum_lycopersicum"

export PATH="$base/bin:$PATH"
export base="$base/$species"
cd $base
if [ "$(pwd)" != "$base" ]; then echo "error, folder not found: $base"; fi
```
Create core folders  
```
cd $base
mkdir -p arx annot data mask
```
### Assembly
Download genomic sequence and reformat it:  
 * simplified FASTA defline with a first word in defline as a unique sequence ID
 * select only nuclear DNA (exclude organelles)
 * set sequence in all uppercase

When possible use genomic sequence from NCBI.  
Match sequence ID in FASTA file with sequence ID in annotation file.  
Use ID from annotation.  
Keep IDs in the file "list.tbl".  
First column in the table is sequence ID and second column is annotation ID.  

Description of assembly is at https://www.ncbi.nlm.nih.gov/assembly/GCF_000188115.4  
Genome assembly version: SL3.0  
* GenBank assembly accession: GCA_000188115.3 (latest)  
* RefSeq assembly accession: GCF_000188115.4 (latest)  
GenBank = Refseq assembly, excluding organelles  
RefSeq annotation ia at https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Solanum_lycopersicum/103/  
Tomato international community project Sol Genomics Network https://solgenomics.net/organism/Solanum_lycopersicum/genome  
```
cd $base/arx
mkdir genbank refseq itag
```
### GenBank assembly
```
cd $base/arx/genbank/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_genomic.fna.gz
gunzip GCA_000188115.3_SL3.0_genomic.fna.gz

cat GCA_000188115.3_SL3.0_genomic.fna | grep '^>' | grep -v AEKE

>CM001064.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 1, whole genome shotgun sequence
>CM001065.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 2, whole genome shotgun sequence
>CM001066.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 3, whole genome shotgun sequence
>CM001067.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 4, whole genome shotgun sequence
>CM001068.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 5, whole genome shotgun sequence
>CM001069.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 6, whole genome shotgun sequence
>CM001070.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 7, whole genome shotgun sequence
>CM001071.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 8, whole genome shotgun sequence
>CM001072.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 9, whole genome shotgun sequence
>CM001073.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 10, whole genome shotgun sequence
>CM001074.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 11, whole genome shotgun sequence
>CM001075.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 12, whole genome shotgun sequence

probuild --stat --seq GCA_000188115.3_SL3.0_genomic.fna --details

LABEL          GCA_000188115.3_SL3.0_genomic.fna
SEQUENCE_SIZE  827747456
SEQUENCE_ACGT  460960521
NT_A           146242095
NT_C           84168461
NT_G           84257462
NT_T           146292503
NT_N           81405486
SEQUENCE_atcg  285381449
NT_a           99612832
NT_c           43075446
NT_g           43110159
NT_t           99583012
NT_n           0
SEQUENCE_other 0
GC             34.1
RECORDS        3148

gzip GCA_000188115.3_SL3.0_genomic.fna
```
### RefSeq 103 assembly
```
cd $base/arx/refseq/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.4_SL3.0/GCF_000188115.4_SL3.0_genomic.fna.gz
gunzip GCF_000188115.4_SL3.0_genomic.fna.gz

grep '^>' GCF_000188115.4
>NC_015438.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 1, SL3.0, whole genome shotgun sequence
>NC_015439.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 2, SL3.0, whole genome shotgun sequence
>NC_015440.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 3, SL3.0, whole genome shotgun sequence
>NC_015441.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 4, SL3.0, whole genome shotgun sequence
>NC_015442.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 5, SL3.0, whole genome shotgun sequence
>NC_015443.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 6, SL3.0, whole genome shotgun sequence
>NC_015444.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 7, SL3.0, whole genome shotgun sequence
>NC_015445.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 8, SL3.0, whole genome shotgun sequence
>NC_015446.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 9, SL3.0, whole genome shotgun sequence
>NC_015447.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 10, SL3.0, whole genome shotgun sequence
>NC_015448.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 11, SL3.0, whole genome shotgun sequence
>NC_015449.3 Solanum lycopersicum cultivar Heinz 1706 chromosome 12, SL3.0, whole genome shotgun sequence
>NC_035963.1 Solanum lycopersicum bio-material TGRC:LA1421 mitochondrion, complete genome
>NC_007898.3 Solanum lycopersicum chloroplast, complete genome

probuild --stat --seq GCF_000188115.4_SL3.0_genomic.fna --details
LABEL          GCF_000188115.4_SL3.0_genomic.fna
SEQUENCE_SIZE  828349174
SEQUENCE_ACGT  461562239
NT_A           146412621
NT_C           84298992
NT_G           84386900
NT_T           146463726
NT_N           81405486
SEQUENCE_atcg  285381449
NT_a           99612832
NT_c           43075446
NT_g           43110159
NT_t           99583012
NT_n           0
SEQUENCE_other 0
GC             34.1
RECORDS        3150
```
### ITAG 3.20 assembly
```
cd $base/arx/itag/
wget ftp://ftp.solgenomics.net/tomato_genome/assembly/build_3.00/S_lycopersicum_chromosomes.3.00.fa

grep '^>' S_lycopersicum_chromosomes.3.00.fa
>SL3.0ch00
>SL3.0ch01 
>SL3.0ch02 
>SL3.0ch03 
>SL3.0ch04 
>SL3.0ch05 
>SL3.0ch06 
>SL3.0ch07 
>SL3.0ch08 
>SL3.0ch09 
>SL3.0ch10 
>SL3.0ch11 
>SL3.0ch12

probuild --stat --seq  S_lycopersicum_chromosomes.3.00.fa  --details
LABEL          S_lycopersicum_chromosomes.3.00.fa
SEQUENCE_SIZE  828076956
SEQUENCE_ACGT  746357470
NT_A           245858945
NT_C           127247685
NT_G           127371487
NT_T           245879353
NT_N           81719486
SEQUENCE_atcg  0
NT_a           0
NT_c           0
NT_g           0
NT_t           0
NT_n           0
SEQUENCE_other 0
GC             34.1
RECORDS        13
```
### Format genomic sequence for project  
Using NCBI RefSeq genome sequence  
Rename chromosome names in ITAG style as: chr_01 .. chr_12.  
```
cd $base/arx
ln -s refseq/GCF_000188115.4_SL3.0_genomic.fna
grep '^>' GCF*.fna  | grep -v '^>NW_' | grep -v NC_035963 | grep -v NC_007898 | cut -f1,8 -d' ' | sed 's/^>//'  | awk '{ printf("%s chr_%02d\n", $1, $2) }' > list_ncbi.tbl

get_fasta_with_tag.pl --swap --in GCF_000188115.4_SL3.0_genomic.fna  --out genome.fasta  --list list_ncbi.tbl --v
probuild --stat --details --seq genome.fasta
probuild --reformat_fasta --in genome.fasta --out ../data/genome.fasta --uppercase 1 --letters_per_line 60 --original
rm genome.fasta
probuild --stat --details --seq ../data/genome.fasta

# NCBI softmasking to GFF coordinates
cd $base/arx/refseq
soft_fasta_to_3 < GCF_000188115.4_SL3.0_genomic.fna | sed 's/ S.*e\t/\t/' | awk '{print $1 "\tNCBI\tRepeat\t" $2+1 "\t" $3 "\t.\t.\t.\trepeat_by_ncbi"  }' > ncbi_softmask.gff

cd $base/arx/refseq
gzip GCF_000188115.4_SL3.0_genomic.fna
cd $base/arx/itag
gzip S_lycopersicum_chromosomes.3.00.fa
```
### Annotation
GenBank GFF3 annotation file was empty. Thus, parsing only NCBI RefSeq and ITAG annoatations.
```
cd $base/arx/refseq
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.4_SL3.0/GCF_000188115.4_SL3.0_genomic.gff.gz
gunzip GCF_000188115.4_SL3.0_genomic.gff.gz

gff_to_gff_subset.pl  --swap  --list ../list_ncbi.tbl  --v  --in GCF_000188115.4_SL3.0_genomic.gff  --out refseq.gff

# move masking coordinates to new seq ID and subset it
gff_to_gff_subset.pl  --swap  --list ../list_ncbi.tbl  --v  --out $base/annot/ncbi_softmask.gff  --in ncbi_softmask.gff

cd $base/arx/itag
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.2_release/ITAG3.2_gene_models.gff
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.2_release/ITAG3.2_REPET_repeats_agressive.gff
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.2_release/ITAG3.2_RepeatModeler_repeats_light.gff

gff_to_gff_subset.pl  --swap  --list ../list_itag.tbl  --v  --in ITAG3.2_gene_models.gff --out itag.gff

# move masking coordinates to new seq ID and subset it
gff_to_gff_subset.pl  --swap  --list ../list_itag.tbl  --v  --out $base/annot/itag_RM_softmask.gff  --in  ITAG3.2_RepeatModeler_repeats_light.gff
gff_to_gff_subset.pl  --swap  --list ../list_itag.tbl  --v  --out $base/annot/itag_REPET_softmask.gff  --in ITAG3.2_REPET_repeats_agressive.gff

gzip ITAG3.2_RepeatModeler_repeats_light.gff
gzip ITAG3.2_REPET_repeats_agressive.gff
```
### Masking
```
cd $base/data

# mask genome using REPET aggressive masking coordinates
/home/tool/bedtools2/bin/bedtools  maskfasta  -fi genome.fasta  -bed ../annot/itag_REPET_softmask.gff  -fo genome_REPET.fasta  -soft

# mask genome using RepeatModeler light masking coordinates
/home/tool/bedtools2/bin/bedtools  maskfasta  -fi genome.fasta  -bed ../annot/itag_RM_softmask.gff  -fo genome_RM.fasta  -soft
```
### RnaSeq
```
cd $base
mkdir rnaseq

cd $base/arx
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR795/SRR7959012/SRR7959012.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR795/SRR7959019/SRR7959019.sra

cd $base/rnaseq

mkdir hisat2db run_hisat2
PATH=$PATH:/storage/braker_c/bin/hisat2/
cd $base/rnaseq/hisat2db

cd $base/rnaseq/run_hisat2
nohup hisat2 -p 15 -x ../genome  -1 $base/arx/SRR7959012_1.fastq -2 $base/arx/SRR7959012_2.fastq -S SRR7959012.sam --novel-splicesite-outfile SRR7959012.introns &
nohup hisat2 -p 15 -x ../genome  -1 $base/arx/SRR7959019_1.fastq -2 $base/arx/SRR7959019_2.fastq -S SRR7959019.sam --novel-splicesite-outfile SRR7959019.introns &

# get intons from hisat2
/storage/braker_c/bin/samtools/samtools view -@10 -bSh SRR7959012.sam -o SRR7959012.bam
/storage/braker_c/bin/samtools/samtools view -@10 -bSh SRR7959019.sam -o SRR7959019.bam

/storage/braker_c/bin/samtools/samtools sort -@10 -n SRR7959012.bam -o SRR7959012.s.bam
/storage/braker_c/bin/samtools/samtools sort -@10 -n SRR7959019.bam -o SRR7959019.s.bam

# parisng of introns is slow !

cd $base/rnaseq
mkdir stardb run_star
cd $base/rnaseq/stardb
/home/tool/STAR/bin/Linux_x86_64/STAR  --runMode genomeGenerate --genomeDir . --genomeFastaFiles $base/data/genome_REPET.fasta  --runThreadN 8
cd $base/rnaseq/run_star
/home/tool/STAR/bin/Linux_x86_64/STAR --genomeDir ../stardb/ --runThreadN 16  --readFilesIn $base/arx/SRR7959012_1.fastq,$base/arx/SRR7959019_1.fastq $base/arx/SRR7959012_2.fastq,$base/arx/SRR7959019_2.fastq

```
