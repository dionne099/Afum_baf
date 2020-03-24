#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Location::Simple;

my $db = 'Afum_genomes.asm.fa';
my $dbh= Bio::DB::Fasta->new($db);

my $infile = shift || "baf-vs-Afum_genomes.FASTA.tab";
open(my $fh => $infile) || die "$infile: $!";

my $out = Bio::SeqIO->new(-format => 'fasta');
my %seen;
while(<$fh>) {
    my ($q,$h,$pid,$identity,$gaps,$gapopen,$qs,$qe,$hitstart,$hitend,$evalue,$bits) = split;
    my ($sp,$scaffold) = split(/\|/,$h);
    if ( ! exists $seen{$sp}->{$q} || $seen{$sp}->{$q}->[2] > $evalue ) {
	$seen{$sp}->{$q} = [$h,Bio::Location::Simple->new(-start => $hitstart, -end => $hitend),$evalue];
    }
}

for my $sp ( keys %seen ) {
    my %done;
    for my $query ( keys %{$seen{$sp}} ) {
	my ($h,$loc, $evalue) = @{$seen{$sp}->{$query}};
	my $skip = 0;
	foreach my $l ( @{$done{$h}} ) {
	    if ( $l->overlaps($loc) ) {
		$l = $l->union($loc);
		$skip = 1;
		last;
	    }
	}
	unless ( $skip ) {
	    push @{$done{$h}}, $loc;
	
#	next if $done{"$h\_$hitstart\_$hitend"}++;
	    my $chunk = $dbh->seq($h,$loc->start => $loc->end);
	    $chunk = Bio::PrimarySeq->new(-id => sprintf("%s_%d_%d",$h,$loc->start,$loc->end), 
					  -seq => $chunk, -description => "hit=$query");
	    $out->write_seq($chunk);
	}
    }
}
