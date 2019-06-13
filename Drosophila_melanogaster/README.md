# Species: _Drosophila_melanogaster_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Environmental variables setup:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Drosophila_melanogaster"

export PATH="$base/bin:$PATH"
export base="$base/$species"
cd $base
if [ "$(pwd)" != "$base" ]; then echo "error, folder not found: $base"; fi
```
Create core folder structure
```
cd $base
mkdir -p arx annot data
```
Download genomic sequence and reformat it:  
 * simplified FASTA defline with unique sequence ID as a first word in defline
 * select only nuclear DNA (exclude organelles)
 * all uppercase
Use genomic sequence from NCBI, when possible.  
```
cd $base/arx

``
