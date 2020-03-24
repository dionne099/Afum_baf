#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;
use Bio::DB::Fasta;

my $db = 'Afum_genomes.asm.fa';
my $dbh= Bio::DB::Fasta->new($db);

my $infile = shift || "baf_mRNA-vs-Afum_genomes.DNA.domtbl";
open(my $fh => $infile) || die "$infile: $!";

my $out = Bio::SeqIO->new(-format => 'fasta');
while(<$fh>) {
	next if /^#/;
		my ($h,undef,$hlen,$q,undef,$qlen,
	    $fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
	    $score,$dombias,
	    $qstart, $qend,$hitstart,$hitend, $envfrom,$envto,$acc,$desc) =
		split(/\s+/,$_,23);
    next if $ivalue > 1e-50;
    my $chunk= $dbh->seq($h,$hitstart => $hitend); 
    $chunk = Bio::PrimarySeq->new(-id => "$h\_$hitstart\_$hitend", -seq => $chunk);
    $out->write_seq($chunk)
}
