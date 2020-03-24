#!/usr/bin/bash
#SBATCH --mem 8gb -p short -n 8 -N 1 --out run_pipeline.log

module load hmmer
module unload perl
module load miniconda3
module load fasta
module load mafft
module load hmmer/3
module load fasttree
fasta36 -E 1e-5 -m 8c baf_mRNA.fa Afum_genomes.asm.fa > baf-vs-Afum_genomes.FASTA.tab
perl extract_DNA_hits_native.pl baf-vs-Afum_genomes.FASTA.tab > baf-vs-Afum_genomes.FASTA.seqs
ssearch36 -E 1e-20 -m 8c baf-vs-Afum_genomes.FASTA.seqs baf_mRNA.fa > bafHits-to-names.SSEARCH.tab
ssearch36 -E 1e-20 baf-vs-Afum_genomes.FASTA.seqs baf_mRNA.fa > bafHits-to-names.SSEARCH
perl rename_hits.pl baf-vs-Afum_genomes.FASTA.seqs bafHits-to-names.SSEARCH.tab > baf-vs-Afum_genomes.FASTA.seqs_renamed

cat baf_mRNA.fa baf-vs-Afum_genomes.FASTA.seqs_renamed | hmmalign baf_mRNA.hmm - >baf.hmmalign.msa
esl-reformat --replace=x:- --gapsym=- afa baf.hmmalign.msa > baf.hmmalign.aln

cat baf_mRNA.fa baf-vs-Afum_genomes.FASTA.seqs_renamed | mafft - > baf.mafft.aln

FastTreeMP -gtr -gamma < baf.mafft.aln > baf.mafft.tre
FastTreeMP -gtr -gamma < baf.hmmalign.aln > baf.hmmalign.tre
