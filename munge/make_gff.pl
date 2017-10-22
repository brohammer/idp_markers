#!/usr/bin/perl -w 
use strict;
use Getopt::Std; 
## About: This script will produce a BED file from the blast output list. This can then be used to extract the target 
# sequence between the mapped markers for pairwise alignment. 

our( $opt_i, $opt_o ); 
getopts("i:o:");

open ( my $infile, '<', $opt_i ); 
open ( my $outfile, '>', $opt_o ); 

# Loop through blast file and retrieve coordinates and other relevant info
while ( my $line = <$infile> ) {
	chomp $line;
	my ( $id, undef, $forward, $forward_loc, $reverse, $reverse_loc ) = split '\t', $line; 

	my ( $chr, $fstart, $fstop, $strand ); 
	
	# Parse mapping location strings got chr, start, stop, and strand info
	if ( $forward_loc =~ /((chr)?\d+):(\d+)-(\d+);(.*)/ ) { 
		$chr = $1; 
		$fstart = $3-1; 
		$fstop = $4; 
                $strand = $5;
	}

	my ( $rstart, $rstop );
	if ( $reverse_loc =~ /(chr)?\d+:(\d+)-(\d+);/ ) { 
		$rstart = $2-1;
		$rstop = $3; 
	}

	# If the target sequence is greater than 10kb then ignore it
	if ( abs($fstart - $rstart) > 10000 ) {
		warn "skipping $id\n";
		next;
	}

	# Make some modifications for unequal strand mappings  
	if ( $strand eq -1 ) {
		$strand = "-";
	}
	else {
		$strand = "+";
	}

	if ( $fstart > $rstart ) {
		print $outfile "$chr\t$rstart\t$fstop\t$id\t.\t$strand\n";
	}
	else {
		print $outfile "$chr\t$fstart\t$rstop\t$id\t.\t$strand\n";
	}
}

close $infile; 
close $outfile; 
exit 0; 

