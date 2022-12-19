#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import ete3
import ast
import hashlib

workdir = 'output_test'
### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('--node_basename',default=None,help='Input a basename of nodes in output relations. Each node will be enumerated using the basename as prefix.\nIf unset, will generate a node name hash based on leaf names (should ensure node unique naming in any modified database)')
parser.add_argument('-r','--root_name',default=None,nargs='+',help='Specify the root name in output relations')
parser.add_argument('--remove_outgroup',help='If specified, will remove the outgroup sample from output relations',action='store_true')
#/
# parse input
#args = parser.parse_args(['-wd','tularensis_realtry1','--remove_outgroup','--root_name','Francisella tularensis','--node_basename','node_'])
args = parser.parse_args()
workdir = args.workdir
node_basename = args.node_basename
root_name = ' '.join(args.root_name)
remove_outgroup = args.remove_outgroup
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (created by metadata-command). Please check the input:')
    print(workdir)
    sys.exit()
if not os.path.exists(workdir+'/'+'genomes_derep_representants.dnd'):
    print('Could not locate dereplicated genome representants tree-file (.dnd)!')
    print('Please make sure you have run the phylo-command on the workdir!')
    sys.exit()
if not os.path.exists(workdir+'/'+'outgroup_accession.txt'):
    print('Could not locate file outgroup_accession.txt to use for tree rooting')
    sys.exit()
#/
###/

### Parse selected accessions metadata
try:
    print('Parsing metadata...')
    with open(workdir+'/'+'metadata_selected.tsv','r') as f:
        line = f.readline()
        accessions = ast.literal_eval(line)
except:
    print('Could not read metadata at location:')
    print(workdir+'/'+'metadata_selected.tsv')
    print('Have you run the metadata-command?')
    sys.exit()
###/

### Parse tree
tree = None
with open(workdir+'/'+'genomes_derep_representants.dnd','r') as f:
    tree = ete3.Tree(f.readline().strip('\n'))
if not tree:
    print('Tree-file was empty or not located')
    sys.exit()
###/

### Root tree by outgroup
outgroup_accession = None
with open(workdir+'/'+'outgroup_accession.txt','r') as f:
    outgroup_accession = f.readline().strip('\n')
###/

### Add node names, get leaves (samples), and get leaves path to root
tree.set_outgroup(outgroup_accession)

# Get leaves and set names to nodes that are not named
node_name_iter = 0
leaves_nodes = {}
for node in tree.iter_descendants():
    if not node.name:
        # Check if we have user-specified basename to nodes
        if node_basename:
            node.name = node_basename+str(node_name_iter)
            node_name_iter += 1
        else: # else, assign a random hash (based on leaves to node) as node name
            node_leaves = []
            for entry in node.get_leaves():
                node_leaves.append(entry.name)
            node.name = hashlib.md5(' '.join(node_leaves).encode('utf-8')).hexdigest()
        #/
    else:
        leaves_nodes[node.name] = []
#/
# Parse paths from each leaf to root
for node in tree.iter_descendants():
    if node.name in leaves_nodes:
        for ancestor in node.get_ancestors():
            # check if ancestor, then modify name
            if ancestor.is_root():
                ancestor.name = 'root'
            #/
            leaves_nodes[node.name].append(ancestor.name)
#/
###/

### Remove outgroup and set root
if remove_outgroup:
    del leaves_nodes[outgroup_accession] # delete outgroup
    
    # set root of leaves
    for leaf,nodes in leaves_nodes.items():
        del nodes[-1] # remove "root"
        
        # rename branch at root
        if root_name:
            nodes[-1] = root_name
        else:
            nodes[-1] = 'root'
        #/
    #/
###/

### Compile output parent-child table from leaves path-to-root and leaves to genome files
## Parent-child file
with open(workdir+'/'+'derep_genomes_tree2tax.tsv','w') as nf:
    # write header
    writeArr = ['child','parent']
    nf.write('\t'.join(writeArr)+'\n')
    #/
    # write leaf parent-child relations from leaf to root
    written_arrs = []
    for leaf,nodes in leaves_nodes.items():
        path = [leaf] + nodes
        for i in range(len(path)):
            if i >= len(path)-1: break # stop when we are out of range
            child = path[i]
            parent = path[i+1]
        
            writeArr = [child,parent]
            if writeArr in written_arrs: continue # skip if redundant row
            nf.write('\t'.join(writeArr)+'\n')
            
            written_arrs.append(writeArr) # keep track of redundant row
    #/
##/
## Leaf-genome file
with open(workdir+'/'+'derep_genomes_map.tsv','w') as nf:
    for leaf in leaves_nodes:
        name_split = leaf.split('_')
        accession = name_split[-2]+'_'+name_split[-1] # GCx_0000000
        writeArr = [accession,leaf]
        nf.write('\t'.join(writeArr)+'\n')
##/
###/

###
"""
IDE-guide for flextaxd (run on Carl):
    
copy original database:
    cp /home/databases/kraken/kraken2/NCBI_GTDB_2022_07_r207_ftul/sources/databases/francisellaceae.db .

visualise original database:
    flextaxd -db francisellaceae.db --vis_type tree --vis_depth 0 --visualise_node Francisellaceae
    
run my pipeline, at output/tree2tax-script, specify --root_name B.137 --remove_outgroup

modify original database using flextaxd:
    flextaxd -db francisellaceae.db --mod_file ../tularensis3/derep_genomes_tree2tax.tsv --genomeid2taxid ../tularensis3/derep_genomes_map.tsv --replace --parent B.137
    
visualize modified database:
    flextaxd -db francisellaceae.db --vis_type tree --vis_depth 0 --visualise_node B.135

visualize full database:
    flextaxd -db francisellaceae.db --vis_type tree --vis_depth 0 --visualise_node Francisellaceae
"""
###/