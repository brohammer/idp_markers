#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Bio::Index::Fasta;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Data::Printer;
##About: This script is used to blast genes that did not have a reciprocal match against the chromosome space

our ($opt_i,$opt_s,$opt_b);
getopts("i:s:b:");

#my $genomeDB = '/home/maize/shared/databases/blast/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa';
my $genomeDB = '/home/maize/shared/databases/blast/Zea_mays/PH207/ZmaysPH207_443_v1.0.fa';
my $indexName = 'idpPrimerSeqs';
my $fastaIndex = Bio::Index::Fasta->new( -filename => $indexName,
										 -write_flag => 1);
$fastaIndex->make_index($opt_s);

my $binRef = &fetchTargets($opt_b);

open ( my $infile, '<', $opt_i );
while ( my $line = <$infile> ) {
  chomp $line;
  my ( $bin, $set, $p1, $seq1, $p2, $seq2 ) = split '\t', $line; 
  my $fastaID1 = join(";", $p1,$bin,$set);
  my $fastaID2 = join(";", $p2,$bin,$set);
  my $results1 = &runBlast($fastaIndex, $fastaID1, $genomeDB);
  my $results2 = &runBlast($fastaIndex, $fastaID2, $genomeDB);
  
  if ( $results1 eq 'FAIL' || $results2 eq 'FAIL' ) {
	next;
  }

  my $parsedResult1 = &parseBlast($$results1);
  my $parsedResult2 = &parseBlast($$results2);

  if  ( $parsedResult1 ne 'FAIL' && $parsedResult2 ne 'FAIL' ) {
		my $toPrint = join("\t", $set,$bin,$p1,$parsedResult1,$p2,$parsedResult2);
		print $toPrint,"\n";
  }
	`rm $$results1`;
	`rm $$results2`;
}
close $infile;
exit 0;

sub fetchTargets {
	open( my $infile, '<', $_[0] );
	my %hash; 
	while ( my $line = <$infile> ) {
		chomp $line; 
		my ( $bin, $start, $stop ) = split '\t', $line; 
		my $chr; 
		if ( $bin =~ /(\d+)\.\d+/ ) {
			$chr = $1; 
		}

		if ( $chr eq '10' ) {
			$chr = "chr" . $chr;
		}
		else {
			$chr = "chr0" . $chr;
		}

		$hash{$bin}{'chr'} = $chr; 
		$hash{$bin}{'start'} = $start-10000;
		$hash{$bin}{'stop'} = $stop+10000; 
	}
	return \%hash;
}
sub runBlast {
    my ( $index, $gene, $db ) = @_;
    my $seq = $index->fetch($gene); # fetch gene sequence from fasta index
	
    if ( !$seq ) { 
	return 'FAIL';
    }
    my $factory = Bio::Tools::Run::StandAloneBlastPlus->new( -db_name => "$db" ) or die "couldn't get blast factory\n";
    # Generate a random string for the temporary output file name
    my $temp_file;
    my @chars = ("A".."Z", "a".."z");
    $temp_file .= $chars[rand @chars] for 1..8;

    # here we call blastn
    $factory->blastn(
        -outfile=> $temp_file,
        -query=> $seq,
        -outformat=> 5,
        -method_args=> [ '-task'=>'blastn', '-num_threads'=>10, '-max_target_seqs'=>1, '-max_hsps'=>1, '-perc_identity'=>100, '-evalue'=>1000, '-ungapped'=>"", '-word_size'=>7, '-reward'=>1, '-penalty'=>-5, '-gapopen'=>3, '-gapextend'=>3]),
        $factory->cleanup;
        return \$temp_file;
}

sub parseBlast {
  my $in = new Bio::SearchIO(-format => 'blastxml', -file => $_[0] );

  # Loop through results
  my $result = $in->next_result;

	# Get hits as long as they are defined
  if ( $result->num_hits < 1 ) {
    return 'FAIL';
  }

  while ( my $hit = $result->next_hit ) {
  	my $hsp = $hit->next_hsp;

		my ( $indMarker, $bin, $markerSet ) = split ';', $result->query_description;

		my $chrom = $binRef->{$bin}->{'chr'};
		my $startBoundary = $binRef->{$bin}->{'start'};
		my $stopBoundary = $binRef->{$bin}->{'stop'};
		
		my $locInfo = $hit->name . ":" . $hsp->hit->start . "-" . $hsp->hit->end . ";" . $hsp->hit->strand;

		my $chr_condition = $hit->name eq $chrom;
		my $start_condition = ($hsp->hit->start < $stopBoundary) && ($hsp->hit->start > $startBoundary);
		my $stop_condition = ($hsp->hit->end < $stopBoundary) && ($hsp->hit->end > $startBoundary);
		
		if ( $chr_condition ) {
			if ( $start_condition || $stop_condition ) {
				return $locInfo;
			}
		}
	}	
  return 'FAIL';
}
