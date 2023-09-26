#!python

import os
import sys
import subprocess
import argparse
import ete3
import ast
import hashlib

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('--node_basename',default=None,help='Input a basename of nodes in output relations. Each node will be enumerated using the basename as prefix.\nIf unset, will generate a node name hash based on leaf names (should ensure node unique naming in any modified database)')
parser.add_argument('-r','--root_name',default=None,nargs='+',help='Specify the root name in output relations')
parser.add_argument('--remove_outgroup',help='If specified, will remove the outgroup sample from output relations',action='store_true')
parser.add_argument('--all_genomes',action='store_true',help='If specified, will run on all genomes and not on de-replicated genomes')
#/
# parse input
args = parser.parse_args()
workdir = args.workdir
node_basename = args.node_basename

root_name = 'root'
if args.root_name:      root_name = ' '.join(args.root_name).replace('"','')

remove_outgroup = args.remove_outgroup
run_on_all_genomes = args.all_genomes
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (created by metadata-command). Please check the input:')
    print(workdir)
    sys.exit()
if not run_on_all_genomes and not os.path.exists(workdir+'/'+'genomes_derep_representants.dnd'):
    print('Could not locate dereplicated genome representants tree-file (.dnd)!')
    print('Please make sure you have run the phylo-command on the workdir!')
    sys.exit()
if run_on_all_genomes and not os.path.exists(workdir+'/'+'genomes.dnd') and not run_on_all_genomes:
    print('Could not locate genomes tree-file (.dnd)!')
    print('Please make sure you have run the phylo-command on the workdir (using the --all_gnomes flag)!')
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

### Set input/output paths depending on input type (all genomes vs. dereplicated genomes)
# Check if run on de-rep genomes
if not run_on_all_genomes:
    tree_input_file = workdir+'/'+'genomes_derep_representants.dnd'
    output_tree2tax_file = workdir+'/'+'derep_genomes_tree2tax.tsv'
    output_genomes_map_file = workdir+'/'+'derep_genomes_map.tsv'
#/
# Else, run on all genomes
else:
    tree_input_file = workdir+'/'+'genomes.dnd'
    output_tree2tax_file = workdir+'/'+'genomes_tree2tax.tsv'
    output_genomes_map_file = workdir+'/'+'genomes_map.tsv'
#/
###/

### Parse tree
tree = None
with open(tree_input_file,'r') as f:
    tree = ete3.Tree(f.readline().strip('\n'))
if not tree:
    print('Tree-file was empty or not located')
    sys.exit()
###/

### Root tree by outgroup
# Get accession
outgroup_accession = None
with open(workdir+'/'+'outgroup_accession.txt','r') as f:
    outgroup_accession = f.readline().strip('\n')
#/
# Find file in outgroup directory
outgroup_file = None
for file_ in os.listdir(workdir+'/'+'outgroup'):
    if file_.find(outgroup_accession) != -1:
        outgroup_file = file_.strip('.fasta')
#/
# Print info
print('Using '+outgroup_file+' as outgroup')
#/
###/

### Add node names, get leaves (samples), and get leaves path to root
tree.set_outgroup(outgroup_file)

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
    del leaves_nodes[outgroup_file] # delete outgroup
    
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
with open(output_tree2tax_file,'w') as nf:
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
with open(output_genomes_map_file,'w') as nf:
    for leaf in leaves_nodes:
        name_split = leaf.split('_')
        accession = name_split[-2]+'_'+name_split[-1] # GCx_0000000
        writeArr = [accession,leaf]
        nf.write('\t'.join(writeArr)+'\n')
##/
###/
