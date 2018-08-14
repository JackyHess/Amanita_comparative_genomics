# Amanita_comparative_genomics

Scripts:

1) Amanita_geneSpaceExpansionPermut.py 

Permutation sampling to determine the expected distribution of copy number (CN) increase on shared sets of branches if duplications and losses are placed among gene families at random.

Required arguments (Inference of duplication, losses and copy numbers according to phylogenomics output by NOTUNG 2.8):

GeneCNs.tab:      inferred CN of each gene family on each branch  
tax_map.txt:      taxonomy information mapping parent and child nodes (see example file)  
Origin.tab:       inferred number of gene families originating on each branch  
Duplications.tab: inferred number of duplications in each gene family on each branch  
Losses.tab:       inferred number of losses in each gene family on each branch  
[Num samples]:      Number of samples to run (integer value)  
GeneCNchange.tab  Inferred CN change of each gene family on each branch  
[output prefix]   Prefix for name of the output file  

2)  Amanita_node_ages.py

Script to estimate the age of each duplication

Required arguments:

Treefile:   Newick format (species label in three letter prefix)  
Outfile:    Name of the output file  

3) Amanita_paralog_clusters.py

Script to determine the extent of syntenic clustering of member of the same gene family along scaffolds

Required arguments:

Genes.bed:  BED files with gene predictions for target species (full gene coordinates only, one per gene)  
MCL_clust:  Clustering output (tab-delimited): GeneName	ClusterName	Species  
Species:    Species description matching species in cluster file^^  
Outfile:    Name of output file  
