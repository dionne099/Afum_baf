#!/usr/bin/env perl
use strict;
use warnings;

use Bio::SeqIO;

my $db = 'Afum_genomes.asm.fa';
if ( ! -f "$db.ssi" ) {
    `esl-sfetch --index $db`;
}
my $infile = shift || "baf-vs-Afum_genomes.FASTA.tab";
open(my $fh => $infile) || die "$infile: $!";

my $out = Bio::SeqIO->new(-format => 'fasta');
while(<$fh>) {
    my ($q,$h,$pid,$identity,$gaps,$gapopen,$qs,$qe,$hitstart,$hitend,$evalue,$bits) = split;
    open(my $ifh => "esl-sfetch $db '$h' |") || die $!;
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $ifh);
    while(my $seq = $seqio->next_seq ) {
	my $chunk;
	if ( $hitstart < $hitend ) {
		$chunk = $seq->subseq($hitstart,$hitend);
	} else {
	    $chunk = $seq->trunc($hitend,$hitstart)->revcom->seq;
	}
	$chunk = Bio::PrimarySeq->new(-id => "$h\_$hitstart\_$hitend", -seq => $chunk, -description => "hit=$q");
	$out->write_seq($chunk)
    }
}
