#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 13:25:31 2018

@author: jacky hess

Permutation test for parallel gene space expansion

"""

import pandas as pd
import numpy as np
import random
import copy
import scipy.stats


from collections import Counter


import datetime



def read_tax_map(map_file):
    
    tax_map = {}
    fp_in = open(map_file, "r")
    lines = fp_in.readlines()
    fp_in.close()
    
    for line in lines:
        child, parent = line.strip().split()
        tax_map[child] = parent
    return tax_map

def create_gene_fam_pool(count_file, spec_name):
    
    family_pool = []
    
    for x in range(0,len(count_file)):
        for y in range(0,count_file[spec_name][x]):
            family_pool.append(count_file['GeneFam'][x])
    return family_pool

def sample(genefam_pool, child, num_dup, num_loss, geneFamInd, origin_fam):
    
    # modify the sample pool in each iteration in order to account for growing 
    # and shrinking gene families
    
    fams_gain = []
    gene_fams = copy.deepcopy(genefam_pool)
    gene_fams.extend(origin_fam)
    for i in range(0,num_dup):
        fams_gain.append(random.sample(gene_fams, 1)[0])
        # attach to sampling pool
        gene_fams.append(fams_gain[-1])
    
    fam_count_gain = Counter(fams_gain)
    
    fams_loss = []
    for i in range(0,num_loss):
        fams_loss.append(random.sample(gene_fams, 1)[0])
        # remove sampling pool
        gene_fams.remove(fams_loss[-1])
    
    fam_count_loss = Counter(fams_loss)
    
    fam_df = pd.DataFrame(geneFamInd)
    
    fam_df[child] = 0
    
    for fam in fam_count_gain.keys():
        fam_df.loc[fam_df['GeneFam'] == fam,child] += fam_count_gain[fam]
    
    for fam in fam_count_loss.keys():
        fam_df.loc[fam_df['GeneFam'] == fam,child] -= fam_count_loss[fam]
    
    return fam_df


def calculate_splits(CNchange_table, splits_dict):
    
    # calculate magnitude of expansion in predetermined splits we're interested in
    # normalize by the number of branches in the split
    
    Abr = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr'].append(sum(Abr['Abr']))
    Amu = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Amu'].append(sum(Amu['Amu']))
    Apo = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Apo'].append(sum(Apo['Apo']))
    AbrApo = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['AbrApo'].append(sum(AbrApo['AbrApo']))
    EM = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['EM'].append(sum(EM['EM']))
    Amu_EM = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Amu_EM'].append((sum(Amu_EM['Amu']) + sum(Amu_EM['EM']))/2)
    Abr_Amu = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_Amu'].append((sum(Abr_Amu['Abr']) + sum(Abr_Amu['Amu']))/2)
    Abr_EM = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Abr_EM'].append((sum(Abr_EM['Abr']) + sum(Abr_EM['EM']))/2)
    Abr_AbrApo = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_AbrApo'].append((sum(Abr_AbrApo['Abr']) + sum(Abr_AbrApo['AbrApo']))/2)
    Abr_Apo = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_Apo'].append((sum(Abr_Apo['Abr']) + sum(Abr_Apo['Apo']))/2)
    Apo_Amu = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Apo_Amu'].append((sum(Apo_Amu['Apo']) + sum(Apo_Amu['Amu']))/2)
    Apo_EM = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Apo_EM'].append((sum(Apo_EM['Apo']) + sum(Apo_EM['EM']))/2)
    Apo_AbrApo = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Apo_AbrApo'].append((sum(Apo_AbrApo['Apo']) + sum(Apo_AbrApo['AbrApo']))/2)
    AbrApo_Amu = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['AbrApo_Amu'].append((sum(AbrApo_Amu['AbrApo']) + sum(AbrApo_Amu['Amu']))/2)
    Abr_Amu_EM = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Abr_Amu_EM'].append((sum(Abr_Amu_EM['Abr']) + sum(Abr_Amu_EM['Amu']) + sum(Abr_Amu_EM['EM']))/3)
    Abr_Amu_Apo = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_Amu_Apo'].append((sum(Abr_Amu_Apo['Abr']) + sum(Abr_Amu_Apo['Amu']) + sum(Abr_Amu_Apo['Apo']))/3)
    Apo_Amu_EM = CNchange_table[(CNchange_table['Abr'] <= 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Apo_Amu_EM'].append((sum(Apo_Amu_EM['Apo']) + sum(Apo_Amu_EM['Amu']) + sum(Apo_Amu_EM['EM']))/3)
    Abr_Apo_AbrApo = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_Apo_AbrApo'].append((sum(Abr_Apo_AbrApo['Abr']) + sum(Abr_Apo_AbrApo['Apo']) + sum(Abr_Apo_AbrApo['AbrApo']))/3)
    Abr_Amu_AbrApo = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] <= 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_Amu_AbrApo'].append((sum(Abr_Amu_AbrApo['Abr']) + sum(Abr_Amu_AbrApo['Amu']) + sum(Abr_Amu_AbrApo['AbrApo']))/3)
    Abr_Apo_EM = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] <= 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Abr_Apo_EM'].append((sum(Abr_Apo_EM['Abr']) + sum(Abr_Apo_EM['Apo']) + sum(Abr_Apo_EM['EM']))/3)
    Abr_Amu_Apo_EM = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] <= 0) & (CNchange_table['EM'] > 0)]
    splits_dict['Abr_Amu_Apo_EM'].append((sum(Abr_Amu_Apo_EM['Abr']) + sum(Abr_Amu_Apo_EM['Amu']) + sum(Abr_Amu_Apo_EM['Apo']) + sum(Abr_Amu_Apo_EM['EM']))/4)
    Abr_Amu_Apo_AbrApo = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] <= 0)]
    splits_dict['Abr_Amu_Apo_AbrApo'].append((sum(Abr_Amu_Apo_AbrApo['Abr']) + sum(Abr_Amu_Apo_AbrApo['Amu']) + sum(Abr_Amu_Apo_AbrApo['Apo']) + sum(Abr_Amu_Apo_AbrApo['AbrApo']))/4)
    all_five = CNchange_table[(CNchange_table['Abr'] > 0) & (CNchange_table['Amu'] > 0) & (CNchange_table['Apo'] > 0) & (CNchange_table['AbrApo'] > 0) & (CNchange_table['EM'] > 0)]
    splits_dict['all_five'].append((sum(all_five['Abr']) + sum(all_five['Amu']) + sum(all_five['Apo']) + sum(all_five['AbrApo']) + sum(all_five['EM']))/5)
    
    return splits_dict


def calculate_splits_pvals(sample_dict, observed_dict, nperm, outfile):
    
    # calculate P-values for overrepresentation
    fp_out = open(outfile, "w")
    for split in observed_dict.keys():
        samples = np.array(sample_dict[split])
        p_val = float(len(samples[samples > observed_dict[split][0]]) +1)/float(nperm + 1)
        fp_out.write(split+'\t'+str(observed_dict[split][0])+'\t'+str(np.mean(samples))+'\t'+str(p_val)+'\n')
    fp_out.close()
    
    splits_df = pd.DataFrame(sample_dict)
    splits_df.to_csv(outfile+".all_data", sep="\t")    
      
    
def run_analysis(counts, map_file, gene_origin, gene_dups, gene_losses, num_permut, CNchange, outfile_pref):
    
    # read gene counts
    count_file = pd.read_table(counts)
    
    # read map file
    tax_map = read_tax_map(map_file)
    
    # read origin
    origin = pd.read_table(gene_origin)
    
    # read duplications
    dups = pd.read_table(gene_dups)
    
    # read losses
    losses = pd.read_table(gene_losses)
    
    # read CN changes
    CNs = pd.read_table(CNchange)
    
    # generate sampling base for each parent 
    genefam_pool = {}
    for parent in tax_map.values():
        if not genefam_pool.has_key(parent):
            genefam_pool[parent] = create_gene_fam_pool(count_file, parent)
    
    splits_dict = {'Abr': [],
                   'Apo': [],
                   'Amu': [],
                   'AbrApo': [],
                   'EM': [],
                   'Amu_EM': [],
                   'Abr_Amu': [],
                   'Abr_EM': [],
                   'Abr_AbrApo': [],
                   'Abr_Apo': [],
                   'Apo_Amu': [],
                   'Apo_EM': [],
                   'Apo_AbrApo': [],
                   'AbrApo_Amu': [],
                   'Abr_Amu_EM': [],
                   'Abr_Amu_Apo': [],
                   'Apo_Amu_EM': [],
                   'Abr_Apo_AbrApo': [],
                   'Abr_Amu_AbrApo': [],
                   'Abr_Apo_EM': [],
                   'Abr_Amu_Apo_EM': [],
                   'Abr_Amu_Apo_AbrApo': [],
                   'all_five': [] }
    
    observed_dict = {'Abr': [],
                   'Apo': [],
                   'Amu': [],
                   'AbrApo': [],
                   'EM': [],
                   'Amu_EM': [],
                   'Abr_Amu': [],
                   'Abr_EM': [],
                   'Abr_AbrApo': [],
                   'Abr_Apo': [],
                   'Apo_Amu': [],
                   'Apo_EM': [],
                   'Apo_AbrApo': [],
                   'AbrApo_Amu': [],
                   'Abr_Amu_EM': [],
                   'Abr_Amu_Apo': [],
                   'Apo_Amu_EM': [],
                   'Abr_Apo_AbrApo': [],
                   'Abr_Amu_AbrApo': [],
                   'Abr_Apo_EM': [],
                   'Abr_Amu_Apo_EM': [],
                   'Abr_Amu_Apo_AbrApo': [],
                   'all_five': [] }
    
    
    # run iterations
    for x in range(0,int(num_permut)):
        if x%10 == 0:
            print(datetime.datetime.now())
            print "Round "+str(x)
        perm_df = pd.DataFrame(count_file['GeneFam'])
        for child in tax_map.keys():
            CNsample = sample(genefam_pool[tax_map[child]], child, sum(dups[child]) , sum(losses[child]), count_file['GeneFam'], origin[origin[child] > 0]['GeneFam'])
            perm_df = pd.merge(perm_df, CNsample, on='GeneFam')
        # update splits count
        splits_dict = calculate_splits(perm_df, splits_dict)
    # calculate P-values for splits
    observed_dict = calculate_splits(CNs, observed_dict)
    calculate_splits_pvals(splits_dict, observed_dict, int(num_permut), outfile_pref+"_bySplit.tab")
    




if __name__ == '__main__':
    
    import sys
    
    if len(sys.argv) == 9:
        run_analysis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
    else:
        print "run_analysis(counts, map_file, origins, gene_dups, gene_losses, num_permut, CNchange, outfile_pref)"
    