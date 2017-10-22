#!/usr/bin/perl -w 
use strict;
use Getopt::Std; 

our( $opt_i, $opt_o ); 
getopts("i:o:");

open ( my $infile, '<', $opt_i ); 
open ( my $outfile, '>', $opt_o ); 

while ( my $line = <$infile> ) {
	chomp $line;
	my ( $id, undef, $forward, $forward_loc, $reverse, $reverse_loc ) = split '\t', $line; 

	my ( $chr, $fstart, $fstop, $strand ); 
	if ( $forward_loc =~ /(chr\d+):(\d+)-(\d+);(.*)/ ) { 
		$chr = $1; 
		$fstart = $2-1; 
		$fstop = $3; 
                $strand = $4;
	}

	my ( $rstart, $rstop );
	if ( $reverse_loc =~ /chr\d+:(\d+)-(\d+);/ ) { 
		$rstart = $1-1;
		$rstop = $2; 
	}

	if ( abs($fstart - $rstart) > 10000 ) {
		warn "skipping $id\n";
		next;
	}

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

