# Tomato
_Solanum lycopersicum_
```
base=/storage4/alexl/Solanum_lycopersicum

cd /storage4/alexl
mkdir Solanum_lycopersicum
cd Solanum_lycopersicum
mkdir arx data annot

git clone git@github.gatech.edu:al68/tomato.git
ln -s $base/tomato/bin  $base/bin
```
### install tools 
```
cd $base/bin
# Get from softmasked genome coordinates of lowercase letters in tabular form
# script from internet https://www.biostars.org/p/134868/
vi soft_fasta_to_3.l
flex -o soft_fasta_to_3.yy.c soft_fasta_to_3.l
gcc -o soft_fasta_to_3  soft_fasta_to_3.yy.c
```
### Assembly and anotation
Genome assembly version: SL3.0  
NCBI https://www.ncbi.nlm.nih.gov/assembly/GCF_000188115.4  
* GenBank assembly accession: GCA_000188115.3 (latest)  
* RefSeq assembly accession: GCF_000188115.4 (latest)  
GenBank ~ Refseq assembly  
Tomato international community project Sol Genomics Network https://solgenomics.net/organism/Solanum_lycopersicum/genome  
```
cd $base/arx
mkdir genbank refseq itag
```
### GenBank
```
cd $base/arx/genbank/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_genomic.fna.gz
gunzip GCA_000188115.3_SL3.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_genomic.gff.gz
gunzip GCA_000188115.3_SL3.0_genomic.gff.gz

# Attention, no gene coordinates in Genbank annotation file

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
```
### RefSeq 103
```
https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Solanum_lycopersicum/103/

cd $base/arx/refseq/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.4_SL3.0/GCF_000188115.4_SL3.0_genomic.fna.gz
gunzip GCF_000188115.4_SL3.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.4_SL3.0/GCF_000188115.4_SL3.0_genomic.gff.gz
gunzip GCF_000188115.4_SL3.0_genomic.gff.gz

grep '^>' GCF_000188115.4_SL3.0_genomic.fna | grep -v '^>NW_'
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
### ITAG 3.20
```
cd $base/arx/itag/
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.2_release/ITAG3.2_gene_models.gff
wget ftp://ftp.solgenomics.net/tomato_genome/assembly/build_3.00/S_lycopersicum_chromosomes.3.00.fa
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.2_release/ITAG3.2_REPET_repeats_agressive.gff
wget ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.2_release/ITAG3.2_RepeatModeler_repeats_light.gff

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
### Prepare genome to run
```
cd $base/annot
mkdir genome_and_repeats
cd $base/annot/genome_and_repeats

# mask genome using REPET aggressive masking coordinates
/home/tool/bedtools2/bin/bedtools maskfasta -fi $base/arx/itag/S_lycopersicum_chromosomes.3.00.fa -bed  $base/arx/itag/ITAG3.2_REPET_repeats_agressive.gff  -fo genome_REPET_all.fasta  -soft

# mask genome using RepeatModeler light masking coordinates
/home/tool/bedtools2/bin/bedtools maskfasta -fi $base/arx/itag/S_lycopersicum_chromosomes.3.00.fa -bed  $base/arx/itag/ITAG3.2_RepeatModeler_repeats_light.gff  -fo genome_RM_all.fasta  -soft

# NCBI softmasking to GFF coordinates
$base/bin/soft_fasta_to_3 < $base/arx/refseq/GCF_000188115.4_SL3.0_genomic.fna > ncbi_soft_mask.3
cat ncbi_soft_mask.3  |  sed 's/ S.*e\t/\t/'  | awk '{print $1 "\tNCBI\tRepeat\t" $2+1 "\t" $3 "\t.\t.\t.\trepeat_by_ncbi"  }' > ncbi_soft_mask_all.gff
rm ncbi_soft_mask.3

# exclude unplaced, MT and PLST contigs/chromosomes
# rename chromosome names as: chr_01 .. chr_12
# manually create list of chromosomes with: old_name new_name

cat  list_ncbi.tbl
cat  list_itag.tbl

$base/bin/gff_to_gff_subset.pl  -in ncbi_soft_mask_all.gff                                  --out ncbi_soft_mask.gff       --list list_ncbi.tbl -v
$base/bin/gff_to_gff_subset.pl  -in $base/arx/itag/ITAG3.2_RepeatModeler_repeats_light.gff  --out itag_rm_soft_mask.gff    --list list_itag.tbl -v
$base/bin/gff_to_gff_subset.pl  -in $base/arx/itag/ITAG3.2_REPET_repeats_agressive.gff      --out itag_aggr_soft_mask.gff  --list list_itag.tbl -v

$base/bin/get_fasta_with_tag.pl --swap --in genome_REPET_all.fasta --out genome_REPET.fasta --list list_itag.tbl --v
$base/bin/get_fasta_with_tag.pl --swap --in genome_RM_all.fasta    --out genome_RM.fasta    --list list_itag.tbl --v
$base/bin/get_fasta_with_tag.pl --swap --in $base/arx/refseq/GCF_000188115.4_SL3.0_genomic.fna  --out genome_NCBI.fasta  --list list_ncbi.tbl --v

probuild --stat --seq genome_NCBI.fasta --details
probuild --stat --seq genome_RM.fasta --details
probuild --stat --seq genome_REPET.fasta --details

mv genome_REPET.fasta $base/data
mv genome_RM.fasta $base/data
mv genome_NCBI.fasta $base/data

cd $base/annot
$base/bin/gff_to_gff_subset.pl  --in ../arx/refseq/GCF_000188115.4_SL3.0_genomic.gff   --out refseq.gff   --list genome_and_repeats/list_ncbi.tbl  --v
$base/bin/gff_to_gff_subset.pl  --in ../arx/genbank/GCA_000188115.3_SL3.0_genomic.gff  --out genbank.gff  --list genome_and_repeats/list_genbank.tbl  --v
$base/bin/gff_to_gff_subset.pl  --in ../arx/itag/ITAG3.2_gene_models.gff               --out itag.gff     --list genome_and_repeats/list_itag.tbl  --v

cat itag.gff | grep -P '\tCDS\t' | cut -f1-7 | sort -k1,1 -k4,4n -k5,5n | awk '{print $1,"z",$3,$4,$5,$6,".","0"}' | tr ' ' '\t'| uniq > itag_CDS_LR.gff
cat refseq.gff | grep -P '\tCDS\t' | cut -f1-7 | sort -k1,1 -k4,4n -k5,5n | awk '{print $1,"z",$3,$4,$5,$6,".","0"}' | tr ' ' '\t'| uniq > refseq_CDS_LR.gff

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


