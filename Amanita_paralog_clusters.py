#!/usr/env python

# Jacky Hess; June 2017

# Parse cluster file and for each orthogroup determine if the set of paralogs for a species is part of a gene cluster
# Gene clusters are defined as more than one member of the family being within 50 kb on the same scaffold

# input data: 
# MCL clustering output (tab-delimited): GeneName	ClusterName	Species
# BED files with gene predictions for each species

import pandas as pd

def read_gene_predictions(BED_file):

    pred = pd.DataFrame(columns=['Chr','Start','End','Gene'])
    fp_in = open(BED_file,"r")
    line = fp_in.readline()

    while line:
        elmts = line.strip().split()
        pred = pred.append(pd.Series([elmts[0],int(elmts[1]),int(elmts[2]),elmts[3]], index=['Chr','Start','End','Gene']), ignore_index=True)
        line = fp_in.readline()
    fp_in.close()
    print pred
    return pred



def read_clusters(cluster_mem_file):
    
    fp_in = open(cluster_mem_file, 'r')
    lines = fp_in.readlines()
    fp_in.close()

    clusters = {}

    for line in lines:
        gene,cluster,spec = line.strip().split()
        if not clusters.has_key(cluster):
            clusters[cluster] = {}
        if not clusters[cluster].has_key(spec):
            clusters[cluster][spec] = []
        clusters[cluster][spec].append(gene)
    
    return clusters

def find_gene_clusters(trimmed_df):

    print trimmed_df
    
    gene_clusters = {}
    sorted_df = trimmed_df.sort_values(by=['Chr','Start'])

    for gene in sorted_df['Gene'].values:
        chrom = sorted_df[(sorted_df.Gene == gene)]['Chr'].values[0]
        start = int(sorted_df.loc[(sorted_df.Gene == gene)]['Start'].values[0])
        end = int(sorted_df.loc[(sorted_df.Gene == gene)]['End'].values[0])
        #chrom,start,end = sorted_df[sorted_df.Gene == gene,['Chr','Start','End']]
        if gene_clusters.has_key(chrom):
            for clust in gene_clusters[chrom]:
                if int(start) <= clust[-1][2]+50000:
                    clust.append([gene,int(start),int(end)])
                    break
                else:
                    gene_clusters[chrom].append([[gene,int(start),int(end)]])
        else:
            gene_clusters[chrom] = []
            gene_clusters[chrom].append([[gene,int(start),int(end)]])

    return gene_clusters    


def catalogue_gene_clusters(BED_file, cluster_file, spec, outfile):

    gene_predictions = read_gene_predictions(BED_file)

    clusters = read_clusters(cluster_file)

    fp_out = open(outfile,"w")
    fp_out.write('OG\tCluster_name\tNum_members\tChr\tStart\tEnd\tMax_Chr_Coord\tMembers\n')

    for OG in clusters.keys():
        if clusters[OG].has_key(spec) and len(clusters[OG][spec]) > 0:
            trimmed_df = gene_predictions[gene_predictions['Gene'].isin(clusters[OG][spec])]
            gene_clusters = find_gene_clusters(trimmed_df)

            cluster_count = 1
            for chrom in gene_clusters.keys():
                for clust in gene_clusters[chrom]:
                    print chrom
                    print gene_predictions[(gene_predictions.Chr == chrom)]['End'] 
                    fp_out.write(OG+"\t"+OG+"_"+str(cluster_count)+"\t"+str(len(clust))+"\t"+chrom+"\t"+str(clust[0][1])+"\t"+str(clust[-1][2])+"\t"+str(max(gene_predictions[(gene_predictions.Chr == chrom)]['End'].values))+"\t")
                    for gene in clust:
                        fp_out.write(gene[0]+";")
                    fp_out.write("\n")
                    cluster_count += 1
    fp_out.close()

if __name__ == '__main__':

    import sys

    if len(sys.argv) == 5:
        catalogue_gene_clusters(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print "catalogue_gene_clusters(BED_file, cluster_file, spec, outfile)"
            
            
        
