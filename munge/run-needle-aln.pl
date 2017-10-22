#!/usr/bin/perl -w 
use strict;
use Getopt::Std;
use Bio::DB::Fasta;
use Bio::Index::Fasta;
use File::Temp; 

## About: This script will perform a pairwise alignment of target sequences.  


our ($opt_i, $opt_q, $opt_s);
getopts("i:q:s:");

# Make query index and subject index
my $queryIdx = Bio::Index::Fasta->new( -filename => 'B73index',
                                       -write_flag => 1);
$queryIdx->make_index($opt_q);

my $subjIdx = Bio::Index::Fasta->new( -filename => 'PH207index',
                                      -write_flag => 1);
$subjIdx->make_index($opt_s);

open ( my $infile, '<', $opt_i);

# Loop through lists of markers and align each  
while ( my $line = <$infile> ) {
	chomp $line;
	
	# Create a temporary file that holds the fasta for the b73 and ph207 sequence
	my $temp1 = new File::Temp ( UNLINK => 1 );
	my $temp2 = new File::Temp ( UNLINK => 1 );
	open ( my $seqfh1, '>', $temp1 );
	open ( my $seqfh2, '>', $temp2 );
	
	# capture the name of each marker to fetch sequence from index
	my @fields = split '\t', $line; 
	my $id = $fields[1];
	my $seq1 = $queryIdx->fetch($id);
	my $seq2 = $subjIdx->fetch($id);

	# Now create the info for the temporary fasta files
	print $seqfh1 ">b73_",$seq1->display_id,"\n"; 
	print $seqfh1 $seq1->seq,"\n"; 
	print $seqfh2 ">ph207_",$seq2->display_id,"\n";
	print $seqfh2 $seq2->seq,"\n";
	my $out = $seq1->display_id . "_" . $seq2->display_id;
	
	# Once fasta files to be aligned are created, we can perform the actual alignment
	`needle $temp1 $temp2 -gapextend 0 -auto -aformat=fasta $out`;
	`java -jar ~/software/jvarkit/dist/msa2vcf.jar < $out >> mismatches.vcf`;
	`rm $out`;
	`needle $temp1 $temp2 -gapextend 0 -auto $out`;
}
close $infile;
