# EukSpecies BRAKER2

Tomas Bruna, Katharina J Hoff, Alexandre Lomsadze, Mario Stanke, Mark Borodovsky

Georgia Institute of Technology, Atlanta, Georgia, USA

University of Greifswald, Greifswald, Germany

Reference: [BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database](https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1)

## Overview

Genome and annotation preparation protocol for eukaryotic species used in the BRAKER2 project.

Only shared components are described in this readme.
Species specific commands are in species README files.

### Installation
```
git clone git@github.com:gatech-genemark/EukSpecies-BRAKER2.git

cd EukSpecies
cd bin

# follow installation instructions in "bin" README file
```
### Setup
Project is set in bash shell.

### Genome sequence
Download genomic sequence and reformat it:
 * Input sequence should be in FASTA format.
 * Unique sequence ID should have no white space symbols in it.
 * Simplify FASTA definition line (defline). First word in defline should be a unique sequence identifier.
 * Select only nuclear DNA sequences from genome (exclude organelles).
 * Set all sequence letters into uppercase.

### Annotation
Use genomic sequence from GenBank when possible. GenBank website and sequence accession IDs are usually more stable than genome project websites. Conversely, annotation is more frequently up-to-date at genomic project locations. Download the most reliable annotation.
* Match sequence ID in FASTA file with sequence ID in annotation file.
* Use ID from annotation.
* Keep information about genome sequence ID and annotation sequence ID in the file "list.tbl".
* First column in the "list.tbl" table is sequence ID and second column is annotation ID.
