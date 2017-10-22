#!/bin/bash
module load bedtools
cd /home/hirschc1/shared/projects/sussmann

#perl make_gff.pl -i blast_output2.txt -o blast_output.gff
#bedtools getfasta -name -bed blast_output.gff -fi /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa -fo blast_output.fasta 

perl make_gff.pl -i ph207_blast_output.txt2 -o ph207_blast_output.gff.test
perl make_gff.pl -i b73_blast_output.txt2 -o b73_blast_output.gff.test
bedtools getfasta -s -name -bed ph207_blast_output.gff.test -fi /home/maize/shared/databases/blast/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa -fo ph207_blast_output.fasta.test
bedtools getfasta -s -name -bed b73_blast_output.gff.test -fi /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa -fo b73_blast_output.fasta.test
