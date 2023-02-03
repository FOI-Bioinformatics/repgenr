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
parser.add_argument('--no_outgroup',action='store_true',help='If specified, will not include the outgroup organism into the tree calculation')
parser.add_argument('--all_genomes',action='store_true',help='If specified, will run on all genomes and not on de-replicated genomes')
#/
# parse input
#args = parser.parse_args(['-wd','tularensis'])
args = parser.parse_args()

workdir = args.workdir
skip_outgroup = args.no_outgroup
run_on_all_genomes = args.all_genomes
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (created by metadata-command). Please check the input:')
    print(workdir)
    sys.exit()
if run_on_all_genomes and not os.path.exists(workdir+'/'+'genomes'):
    print('Could not locate genomes!')
    print('Please check the input working directory and confirm that you have run the metadata and genome modules')
    sys.exit()
if not run_on_all_genomes and not os.path.exists(workdir+'/'+'genomes_derep_representants'):
    print('Could not locate dereplicated genome representants!')
    print('Please check the input working directory and confirm that you have run the derep-command')
    sys.exit()
#/
###/

### Set input/output paths depending on input type (all genomes vs. dereplicated genomes)
# Check if run on de-rep genomes
if not run_on_all_genomes:
    genome_files_dir = workdir+'/'+'genomes_derep_representants'
    output_file_path = workdir+'/'+'genomes_derep_representants.dnd'
#/
# Else, run on all genomes
else:
    genome_files_dir = workdir+'/'+'genomes'
    output_file_path = workdir+'/'+'genomes.dnd'
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
    cmd_cp = ['cp',workdir+'/outgroup/'+outgroup_file,genome_files_dir]
    subprocess.call(' '.join(cmd_cp),shell=True)
    #/
###/

### Run tree algorithm
cmd_mashtree = ['mashtree',genome_files_dir+'/*.fasta','>',output_file_path]
subprocess.call(' '.join(map(str,cmd_mashtree)),shell=True)
###/

### Remove outgroup sample
if not skip_outgroup:
    os.remove(genome_files_dir+'/'+outgroup_file)
###/