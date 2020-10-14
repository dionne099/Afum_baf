library(ape)
library(phytools)
library(geiger)
library(tidyverse)
library(RColorBrewer)
library(ggtree)

geneCopies <- read.table("strain_gene_counts_update.tab",header=TRUE,sep="\t",row.names=1) 

#straintree_all <- read.tree("Baf_Strains3.SNP.fasttree.tre")
straintree_all <- read.tree("Baf_Strains3.IQ_rooted.tre")
tipnames <- straintree_all$tip.label

to_drop <- setdiff(straintree_all$tip.label, rownames(geneCopies))
to_drop <- append(to_drop, c("F20140","F12636","F16216","F15390","F21572","F12041","F20451"))
straintree <- drop.tip(straintree_all,to_drop)
tipnames <- straintree$tip.label

difftable <- setdiff(rownames(geneCopies),straintree$tip.label)
difftable
geneCopiesFilter <- filter(geneCopies,rownames(geneCopies) %in% straintree$tip.label)
geneCopiesFilter$bafB <- factor(geneCopiesFilter$bafB)
geneCopiesFilter$bafC <- factor(geneCopiesFilter$bafC)
geneCopiesFilter$bafD <- factor(geneCopiesFilter$bafD)
geneCopiesFilter$bafPseudo <- factor(geneCopiesFilter$bafPseudo)
geneCopiesCount  = as.matrix(geneCopiesFilter)

#dim(geneCopiesCount)

#plot(straintree,use.edge.length = FALSE)
#blues <- colorRampPalette(brewer.pal(9, "Blues"))(4)
#phylo.heatmap(straintree,geneCopiesCount,fsize=c(0.9,0.8,0.8),standardize=TRUE,
#              split=c(2,0.4),
#              colors=blues)

#intree <- read.tree("Baf_Strains3.SNP.fasttree.tre")
#intree <- read.newick("Baf_Strains.IQ_rooted.tre")
#straintree <- drop.tip(intree, to_drop)

p = ggtree(straintree,layout='circular') + 
  geom_tiplab(size=2, align=TRUE, linesize=0.1) 

pcirc <- gheatmap(p, geneCopiesCount, width=1.5, offset=3, legend_title="baf Copies") + 
  scale_fill_brewer(direction = 1) +
  scale_y_continuous(expand=c(0, 7)) + scale_x_ggtree()  + theme_tree2() 
pcirc

ggsave("baf_copies_circle.pdf",pcirc, height=10)

p = ggtree(straintree) + 
  geom_tiplab(size=2, align=TRUE, linesize=0.1)

phmap<-gheatmap(p, geneCopiesCount, offset=0.5, width=1, color="white",
         colnames=FALSE, legend_title="baf Copies") +
  scale_x_ggtree() + scale_y_discrete(expand=c(0, 4)) + 
  theme_tree2() + scale_fill_brewer(direction = 1)
phmap
ggsave("baf_copies.pdf",phmap, height=9)
