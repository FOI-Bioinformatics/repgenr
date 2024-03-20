#!/usr/bin/env python3

import argparse
import os
import sys
import shutil


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = 
'''
Unpack clusters, creating one subdirectory per cluster
''')

parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('--no_representant',required=False,default=False,action='store_true',help='If specified, will not output the representative genome for cluster subdirectories (default: output cluster-representant)')
#/
# parse input
args = parser.parse_args()

workdir = args.workdir
no_representant = args.no_representant
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory at:')
    print(workdir)
    sys.exit()

if not os.path.exists(workdir+'/'+'derep_genomes_summary2.tsv'):
    print('Could not locate file derep_genomes_summary2.tsv inside workdir!')
    print('Ensure that you have run the "derep" module successfully')
    print('Termiating!')
    sys.exit()
#/
###/

### Import summary-file (columns: (1) representative genome, (2) contained genome, (3) "contained genome is representative genome [BOOL]" )
print('Reading summary-file of cluster representants and contained genomes',flush=True)
clusters = {} # representative-genome -> [contained genomes + representative genome]
representants_found = set()
with open(workdir+'/'+'derep_genomes_summary2.tsv','r') as f:
    for enum,line in enumerate(f):
        # skip header line
        if enum == 0: continue
        #/
        # parse line
        line = line.strip('\n')
        line = line.split('\t')
        representative_genome,contained_genome,is_self = line
        #/
        # save all representants to keep track of empty outputs (when using --no_representant)
        representants_found.add(representative_genome)
        #/
        # check if user wants to skip cluster representants
        if no_representant:
            if representative_genome == contained_genome: continue
        #/
        # save
        if not representative_genome in clusters:       clusters[representative_genome] = []
        clusters[representative_genome].append(contained_genome)
        #/

if not len(clusters) == len(representants_found):
    print('Some clusters only had a representative genome. These clusters will not be output',flush=True)
    print('Number of clusters with representant only: '+str(len(representants_found.difference(clusters))),flush=True)
###/

### Perform unpack
unpack_dir = workdir+'/'+'derep_unpacker'
# remove previous unpack directory if it existed
if os.path.exists(unpack_dir):
    print('Previous directory of unpacked genome clusters found. Will remove this directory before proceeding',flush=True)
    shutil.rmtree(unpack_dir)
#/
# init unpack dir
print('Begin unpacking',flush=True)
os.makedirs(unpack_dir)
#/
# perform unpacking
for cluster_rep,cluster_containeds in clusters.items():
    cluster_unpack_dir = unpack_dir+'/'+cluster_rep
    os.makedirs(cluster_unpack_dir)
    for contained_genome in cluster_containeds:
        source_path = workdir+'/'+'genomes'+'/'+contained_genome
        target_path = cluster_unpack_dir+'/'+contained_genome
        
        # check if need to add cwd to get absolute path in case of input relative paths
        if source_path[0] != '/':       source_path = os.getcwd()+'/'+source_path
        if target_path[0] != '/':       target_path = os.getcwd()+'/'+target_path
        #/
        
        shutil.copy2(source_path,target_path)
#/
print('Unpacking finished!',flush=True)
###/
