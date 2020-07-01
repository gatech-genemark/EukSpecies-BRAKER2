# Evaluation of Complete Genes

Many annotations contain a large number of incomplete genes: Genes without a full initial or terminal exon. We want to identify such genes for the purposes of gene prediction algorithms evaluation.

Initially, we were detecting the incomplete genes directly from annotations, based on specific annotation format. For example, ENSEMBL based annotations have a missing start/stop codon for incomplete genes in the `.gtf` annotation, NCBI annotations have a flag `partial=TRUE` in the 9th gff column. 

The obvious drawback of detecting incomplete genes this way is that for each new annotation format, we need a new script. Moreover, some annotations, such as tomato, do not have any information about the status of completeness of a genes in the `.gtf/.gff` annotation files.

Therefore, we created an alternative way of evaluating the Gene completeness. The script `bin/findPartialGenes.py` categorizes genes and transcripts into complete and incomplete groups. A complete transcript is any transcript starting with `ATG` start codon and ending with `TAA`, `TAG` or `TGA`.

Below is an evaluation of how groups complete transcripts identified by findPartialGenes.py and previous methods differ, this documents investigates the differences.

## ENSEMBL based annotations

The difference between ENSEMBL specific script and the general script is very small. Further, the general script seems to actually yield better results.

### D. rerio

* Identified as complete only by ENSEMBL specific script: A small number of genes starting with TTG or CTG.
* Identified a as complete only by the general script: A small number of ATGs which do not have a corresponding start codon in ENSEMBL `.gtf`.

### T. nigroviridis

* Identified as complete only by ENSEMBL specific script: Genes starting with TTG or CTG.
* Identified a as complete only by the general script: Genes with a source `genoscope`. All these genes seem to be missing stop codon in ENSEMBL annotation, thus they were identified as incomplete before.

### R. prolixus
* Identified as complete only by ENSEMBL specific script: Genes starting with TTG or CTG.
* Identified a as complete only by the general script: 19, some missing starts/stops in the `.gtf` annotation.


## NCBI based annotations

The difference between NCBI specific and general script is quite large, especially for spider. Many partial genes in NCBI are partial because the sequence was modified to match an expected genes, all such genes have a partial=true flag and a note similar to:

```
Note=The sequence of the model RefSeq transcript was modified relative to this genomic sequence to represent the inferred CDS: added 68 bases not found in genome assembly;
```

Such modification are not detected by checking for a valid start and stop codon. Therefore, in case of NCBI annotations, it is better to rely on the partial=TRUE flag.

### B. terrestris
* Identified as complete only by NCBI specific script: Just one
* Identified a as complete only by the general script: 196

### X. tropicalis
* Identified as complete only by NCBI specific script: 21 genes with random start codons.
* Identified a as complete only by the general script: 274, out of 537 incomplete.

### P. tepidariorum

* Identified as complete only by NCBI specific script: 0
* Identified a as complete only by the general script: 2011


## Other annotations

For other annotations, there is no script for evaluation of incomplete genes. Below is just a summary of the number of incomplete transcripts found.

### D. melanogaster 

* A small number of non-canonical starts and stops (134 transcripts)
    * 7 are TTG starts
    * 45 are CTG starts

### A. thaliana

* Only 100 incomplete transcripts

### C. elegans: 

* 68 total incomplete
    * 18 CTG starts
    * 20 TTG starts
    * 25 GTG starts

### M. truncantula:

* Only 5 incomplete transcripts

### P. trichocarpa:

* 131 incomplete transcripts, with random codons (6 TTG, 1 CTG start)

### S. lycopersicum:

* 4859 incomplete transcripts
* No non-canonical start codon is overrepresented. No stop codon either.


