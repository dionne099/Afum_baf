#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my ($seqfile,$hits) = @ARGV;

my %hit2name;
open(my $fh => $hits ) || die $!;
while(<$fh>) {
	next if (/^#/);
	my ($q,$h,$pid,@rest) = split;
	next if $rest[-2] > 1e-20;
	unless ( exists $hit2name{$q} ) {
		$hit2name{$q} = $h;
	}
}
my $in = Bio::SeqIO->new(-format => 'fasta', -file => $seqfile);
my $out = Bio::SeqIO->new(-format => 'fasta');
while(my $seq = $in->next_seq) {
	my $id = $seq->display_id;
	if ( my $n = $hit2name{$id} ) {
		$seq->display_id("$n|$id");
	}
	$out->write_seq($seq);
}
