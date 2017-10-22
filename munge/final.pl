#!/usr/bin/perl -w
use strict; 
use Getopt::Std; 

our ($opt_b, $opt_p, $opt_m); 
getopts("b:p:m:");

my $b73ref = &makeKey($opt_b); 
my $ph207ref = &makeKey($opt_p);
my $markerRef = &getMismatches($opt_m); 


## BG842342        1.01  1:7885682-7885701;1   1:7885919-7885938;-1 1:7885682-7885701;1   1:7885919-7885938;-1
#chr1 REF ALT BG842342
#chr1 REF ALT BG842342
## BG842342        1.01  1:7885682-7885701;1   1:7885919-7885938;-1 


# Now loop through mismatches and print out all relevant information
foreach my $marker ( keys %{$markerRef} ) { 
	my $toPrint = $marker . "\t" . $b73ref->{$marker}->{'bin'} . "\t" . $b73ref->{$marker}->{'fpos'} . "\t" . $b73ref->{$marker}->{'rpos'} . "\t" . $ph207ref->{$marker}->{'fpos'} . "\t" . $ph207ref->{$marker}->{'rpos'};
	print "#$toPrint\n";
	foreach my $pos (  sort {$a <=> $b} keys %{ $markerRef->{$marker} } ) { 
		print "$marker\t$pos\t$markerRef->{$marker}->{$pos}->{'ref'}\t$markerRef->{$marker}->{$pos}->{'alt'}\n";
	}
}

# This subroutine will collect all of the mismatches from the temporary vcfs
# produced from parsing needle MSA outputs
sub getMismatches {
	my ( $handle ) = $_[0];
	open ( my $mismatch_file, '<', $handle ); 
	my $id; 
	my %variants; 
	while ( my $line = <$mismatch_file> ) {
		chomp $line; 

		if ( $line =~ /^#CHROM/ ) { 
			my @fields = split '\t', $line; 
			# Get id from vcf
			if ( $fields[9] =~ /b73_(.*)/ ) {
				$id = $1; 
			}
		}
		elsif ( $line =~ /^#/ ) { 
			# skip the headers
			next;
		}
		else {
			# store all mismatches 
			my @fields = split '\t', $line;
			my ( $pos, $ref, $alt ) = @fields[1,3,4]; 
			$variants{$id}{$pos}{'ref'} = $ref;
			$variants{$id}{$pos}{'alt'} = $alt; 
		}
	}
	close $mismatch_file;
	return \%variants; 
}


sub makeKey {

	my ( $handle ) = $_[0]; 
	my %hash;
	open ( my $blastfile, '<', $handle ); 
	while ( my $line = <$blastfile> ) {
		chomp $line; 
		my @fields = split '\t', $line; 
		my ( $key, $bin, $fpos, $rpos ) = @fields[0,1,3,5];
		$hash{$key}{'bin'} = $bin; 
		$hash{$key}{'fpos'} = $fpos;
		$hash{$key}{'rpos'} = $rpos;
	}
	close $blastfile;
	return \%hash; 
}
