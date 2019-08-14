# Scripts for processing eukaryotic species
Alex Lomsadze  
Georgia Institute of Technology  
2019  


### GeneMark-ES/ET
```
cd $base/bin

```
### install third-party tools
```
cd $base/bin
# Get from softmasked genome coordinates of lowercase letters in tabular form
# script from internet https://www.biostars.org/p/134868/
vi soft_fasta_to_3.l
flex -o soft_fasta_to_3.yy.c soft_fasta_to_3.l
gcc -o soft_fasta_to_3  soft_fasta_to_3.yy.c
```
