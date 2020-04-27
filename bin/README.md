# Scripts for processing eukaryotic species
Alex Lomsadze, Tomas Bruna
Georgia Institute of Technology
2020

## Installation Instructions

## GaTech Dependencies

### probuild
Download `probuild` from http://topaz.gatech.edu/GeneMark/license_download.cgi -- select **GeneMark-ES/ET/EP** option. `probuild` is located in the root folder, copy it here, the `bin` folder of this project or add it to your `$PATH`.

## Third-party Tools

The following tools are assumed to be in the `$PATH` variable

* [GenomeTools](http://genometools.org/index.html)
* [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
* [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/)
* [RepeatMasker](http://www.repeatmasker.org/RMDownload.html)

### Coordinates from soft-masked sequence

Script which reports coordinates of lowercase letters in tabular format from softmasked genome in FASTA format. Program from https://www.biostars.org/p/134868/
```
cd path_to_project_location
cd bin
mkdir –p src
cd src

cp ../soft_fasta_to_3.l   .
flex -o soft_fasta_to_3.yy.c soft_fasta_to_3.l
gcc -o soft_fasta_to_3  soft_fasta_to_3.yy.c
mv soft_fasta_to_3  ..
rm soft_fasta_to_3.l soft_fasta_to_3.yy.c
```
