#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import gzip
import ast

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('--no_outgroup',default=False,help='If specified, will not include the outgroup organism into the tree calculation')
#/
# parse input
#args = parser.parse_args(['-wd','tularensis'])
args = parser.parse_args()

workdir = args.workdir
skip_outgroup = args.no_outgroup
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (created by metadata-command). Please check the input:')
    print(workdir)
    sys.exit()
if not os.path.exists(workdir+'/'+'genomes_derep_representants'):
    print('Could not locate dereplicated genome representants!')
    print('Please check the input working directory and confirm that you have run the derep-command')
    sys.exit()
#/
###/

### Copy-in outgroup sample
if not skip_outgroup:
    # Get accession
    outgroup_accession = None
    with open(workdir+'/'+'outgroup_accession.txt','r') as f:
        outgroup_accession = f.readline().strip('\n')
    #/
    # Find file in outgroup directory
    outgroup_file = None
    for file_ in os.listdir(workdir+'/'+'outgroup'):
        if file_.find(outgroup_accession) != -1:
            outgroup_file = file_
    #/
    # Print info
    print('Using '+outgroup_file+' as outgroup')
    #/
    # Do the copy-in
    cmd_cp = ['cp',workdir+'/outgroup/'+outgroup_file,workdir+'/genomes_derep_representants/']
    subprocess.call(' '.join(cmd_cp),shell=True)
    #/
###/

### Run tree algorithm
cmd_mashtree = ['mashtree',workdir+'/genomes_derep_representants/*.fasta','>',workdir+'/'+'genomes_derep_representants.dnd']
subprocess.call(' '.join(map(str,cmd_mashtree)),shell=True)
###/

### Remove outgroup sample
if not skip_outgroup:
    os.remove(workdir+'/genomes_derep_representants/'+outgroup_file)
###/