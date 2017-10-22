# IDP Markers

**About:** Map indel polymorphism markers developed between B73 and Mo17 to PH207. Markers with indels between B73 and PH207 may be used as diagnostic fingerprints between B73 and PH207. 

## Background

Marker information was obtained from excel file. Only markers that were not 'close' or 'too small' were used. 

Bin coordinates were obtained from [maizeGDB Binviewer page](https://maizegdb.org/bin_viewer).  

A parsed version of the excel file was constructed: `data/idp_markers.txt`. 

This file was converted to a fasta file: `data/idp_markers.fasta`.  

## Mapping Markers

Forward and reverse markers were aligned to both B73 and PH207 using blastn. 

Blastn parameters were modified to maximize sensitivity: 
* task: blastn
* max_target_seqs: 1
* max_hsps: 1  
* perc_identity: 100
* evalue: 1000  
* ungapped:  
* word_size: 7  
* reward: 1  
* penalty: -5  
* gapopen: 3  
* gapextend 3  

This script, `munge/blast.pl`, was submitted through a PBS wrapper, `munge/blast.job`.  

Blast outputs: `data/b73_blast_output.txt` and `data/ph207_blast_output.txt`.
