#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess
import argparse
import gzip
import ast

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of total threads to use (default: 16)')
parser.add_argument('--no_outgroup',action='store_true',help='If specified, will not include the outgroup organism into the tree calculation')
parser.add_argument('--all_genomes',action='store_true',help='If specified, will run on all genomes and not on de-replicated genomes')
#/
# parse input
#args = parser.parse_args(['-wd','snippy_tularensis3'])
args = parser.parse_args()
num_threads = args.threads

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

### Get absolute path for workdir, needed for software calls
exec_cwd = os.getcwd()
workdir = exec_cwd + '/' + workdir # get absolute path to workdir
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

### Init phylo workdir (snippy and iqtree produce output files)
phylo_wd = workdir+'/'+'phylo_workdir'
# Remove previous phylo wd if it exists
if os.path.exists(phylo_wd):
    shutil.rmtree(phylo_wd)
#/
# init phy wd
os.makedirs(phylo_wd)
#/
###/

### Setup Snippy input files (it expects genomes as sample_name\tpath_to_file)
snippy_ref = None
with open(phylo_wd+'/'+'input.list','w') as nf:
    # write path to genomes, set first genome as reference
    for enum,file_ in enumerate(sorted(os.listdir(genome_files_dir))):
        # make first file the reference
        if enum == 0:
            snippy_ref = file_
            continue
        #/
        nf.write(file_+'\t'+genome_files_dir+'/'+file_+'\n')
    #/
    # check if add-in outgroup
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
        # Add outgroup to snippy input list
        nf.write(outgroup_file+'\t'+workdir+'/outgroup/'+outgroup_file)
        #/
    #/
###/

### Run Snippy (inside phylo_wd, snippy produces some files)
# Check so snippy_ref was obtained
if not snippy_ref:
    print('Was unable to correctly structure data for Snippy! Was unable to set a reference for snippy. Make sure your input is valid, and if so, check code')
    sys.exit('Terminating..')
#/
## change to phylo_wd and run snippy
os.chdir(phylo_wd)

# make snippy-multi file
cmd_snippy_setup = ['snippy-multi','input.list','--ref',genome_files_dir+'/'+snippy_ref,'--cpus',num_threads,'>','snippy-multi.sh']
subprocess.call(' '.join(map(str,cmd_snippy_setup)),shell=True)
#/
# execute snippy-multi file
cmd_snippy_execute = ['bash','./snippy-multi.sh']
subprocess.call(' '.join(map(str,cmd_snippy_execute)),shell=True)
#/
# clean up "large" files
print('Removing temporary alignment-files')
for file_ in os.listdir(genome_files_dir) + os.listdir(workdir+'/'+'outgroup'):
    if os.path.exists(phylo_wd+'/'+file_):
        shutil.rmtree(phylo_wd+'/'+file_)
#/
###/

### Run IQTREE (inside phylo_wd, iqtree produces some file)
cmd_iqtree = ['iqtree','-T','auto','--threads-max',num_threads,'-o',outgroup_file,'-s','core.full.aln']
subprocess.call(' '.join(map(str,cmd_iqtree)),shell=True)
#iqtree -T auto --threads-max 40 -s snippy_workdir/core.full.aln -o Francisellaceae_Francisella_sp000764555_GCF_000764555.1.fasta
###/

### Save which sample was used as reference in Snippy and write IQTREE's *.treefile to upper workdir with "Reference" replaced to its sample name
with open(workdir+'/'+'tree_reference.txt','w') as nf:
    nf.write(snippy_ref+'\n')

with open(phylo_wd+'/'+'core.full.aln.treefile','r') as f:
    with open(output_file_path,'w') as nf:
        for line in f:
            nf.write(line.replace('Reference:',snippy_ref+':'))
###/
