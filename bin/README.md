# Scripts for processing eukaryotic species
Alex Lomsadze  
Georgia Institute of Technology  
2019  

## Install 
GaTech dependences
* probuild
```
cd path_to_project_location
cd bin
mkdir –p src
cd src
git clone git@github.gatech.edu:al68/probuild.git
cd probuild/src
make
cp probuild ../bin
cd ..
rm –rf probuild
```
## Install third-party tools

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

