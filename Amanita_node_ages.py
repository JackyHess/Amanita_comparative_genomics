#!/usr/env python

# Jacky Hess; June 2017

# For each gene family members of a target species, determine the age of all paralogs 

import dendropy

node_ages_L1 = { 'Abr' : set(['Abr']),
              'Amu' : set(['Amu']),
              'Apo' : set(['Apo']),
              'Ath' : set(['Ath']),
              'Ain' : set(['Ain']),
              'Vvo' : set(['Vvo'])}

node_ages_L2 = {'AbrApo' : set(['Abr','Apo']),
                'AS' : set(['Ath','Ain'])
                }

node_ages_L3 = {'EM' : set(['Abr','Apo','Amu'])}



def is_duplication(node):

    children = node.child_nodes()

    taxa_left = []
    for leaf in children[0].leaf_nodes():
        taxa_left.append(leaf.taxon.label[:3])
    taxa_right = []
    for leaf in children[1].leaf_nodes():
        taxa_right.append(leaf.taxon.label[:3])    
    if set(taxa_left).intersection(set(taxa_right)) == set([]):
        return False
    else:
        return True


def node_age(node):

    leaves = node.leaf_nodes()
    leaf_taxa = []

    for leaf in leaves:
        leaf_taxa.append(leaf.taxon.label[:3])
    
    #find the most shallow MRCA for a set of leaf taxa
    age = ''
    for name in node_ages_L1.keys():
        if set(leaf_taxa).union(node_ages_L1[name]) == node_ages_L1[name]:
            age = name
            return age
    for name in node_ages_L2.keys():
        if set(leaf_taxa).union(node_ages_L2[name]) == node_ages_L2[name]:
            age = name
            return age
    for name in node_ages_L3.keys():
        if set(leaf_taxa).union(node_ages_L3[name]) == node_ages_L3[name]:
            age = name
            return age
    return 'root'

def get_node_ages(node):

    if is_duplication(node):
        age = node_age(node)
        return age
    # unless node is root of the tree, move backwards until duplication is found
    elif node.level() != 0:
        return get_node_ages(node.parent_node)
    # if we arrive at the root node report the age of the orthogroup
    else:
        age = node_age(node)
        return age
        
    

def get_paralog_node_ages(tree_file, outfile):

    tree = dendropy.Tree.get(path=tree_file, schema="newick")

    leaves = tree.leaf_nodes()

    fp_out = open(outfile, "w")
    
    for leaf in leaves:
        spec = leaf.taxon.label[:3]
        age = get_node_ages(leaf.parent_node)
        fp_out.write(tree_file.split("/")[-1]+"\t"+leaf.taxon.label+"\t"+spec+"\t"+age+"\n")
    fp_out.close()


if __name__ == '__main__':

    import sys

    if len(sys.argv) == 3:
        get_paralog_node_ages(sys.argv[1], sys.argv[2])
    else:
        print "get_paralog_node_ages(tree_file, outfile)"
        
