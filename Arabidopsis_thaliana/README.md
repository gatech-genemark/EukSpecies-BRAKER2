# Species: _Arabidopsis thaliana_
Alex Lomsadze  
Georgia Institute of Technology  
2019  
## Project setup
Project is set in bash shell.  

Environmental variables setup:  
```
umask 002

base="/storage3/w/alexl/EukSpecies"
species="Arabidopsis_thaliana"

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


