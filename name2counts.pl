#!/usr/bin/env perl
use Bio::SeqIO;
use strict;
use warnings;
my $in = "baf-vs-Afum_genomes.FASTA.seqs_renamed";

my %counts;
my %genes;
my $seqio = Bio::SeqIO->new(-format => 'fasta', -file => $in);
while(my $seq = $seqio->next_seq ){
    my ($id) = $seq->display_id;
    my @idrow = split(/\|/,$id);
    my ($g,$strain,$loc) = @idrow;
    $genes{$g}++;
    $counts{$strain}->{$g}++;
}
my @genes = sort keys %genes;
print join("\t", qw(taxa),@genes),"\n";
for my $strain ( sort keys %counts ) {
    print join("\t", $strain, map { $counts{$strain}->{$_} || 0 } @genes),"\n";
}
