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

## Convert coordinates  

Use the blast hit coordinates to get the target sequence between the forwared and reverse primers:  
* `munge/make_bed.pl` - This script will convert the blast output to a bed file .
*  `munge/get_fasta.pl` - Use this script to call bedtools and get fasta from BED  

The previous scripts can be ran together in the `munge/get_fasta.sh` wrapper. 

## Pairwise alignment  

With both target sequences in hand, now we need to align the B73 and PH207 sequences to see which markers have indels and might make good diagnostic markers.  

I am using the Emboss needle program to perform the alignment:   
* `munge/run-needle-aln.pl`  

This script makes a system call to Java program called `msa2vcf`, which is a part of the jvarkit (~/software/jvarkit/dist/msa2vcf.jar). This program will parse the output of the MSA output and produce a VCF of variant sites. 

## Parsed output  

This vcf output is pretty ugly and hard to read. Use the `munge/pretty_mismatch.pl` script to make this output more readable. 

### Program versions  
* `Emboss 6.5.7`  
* `Blast 2.6.0`  
* `Bioperl 1.7.1`  
* `jvarkit`  
* `BEDtools 2.25.0`  



