#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-t','--threads',type=int,default=24,help='Number of threads to use')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
#/
# parse input
args = parser.parse_args()

workdir = args.workdir
num_threads = args.threads
keep_files = args.keep_files
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (make sure you have run metadata and genome commands):')
    print(workdir)
    sys.exit()
#/
###/

### Check so genomes-dir exists (required by this script)
if not os.path.exists(workdir+'/'+'genomes'):
    print('Could not locate genome directory inside specified workdir (make sure you have run "metadata" and "genome" commands):')
    print(workdir)
    sys.exit()
###/

## Run dRep compare. It will produce a PDF-output with a mash-tree dendogram.
print('Computing ANI mash-tree with dRep using subcommand compare...')
cmd_drep_compare = ['dRep','compare',
                    '--SkipSecondary',
                    '-g',workdir+'/genomes/*.fasta','--processors',num_threads,
                    workdir+'/'+'glance_wd']
subprocess.call(' '.join(map(str,cmd_drep_compare)),shell=True)

if not os.path.exists(workdir+'/glance_wd'):
    print('Could not locate "dRep compare" directory! Please confirm that dRep is installed and can execute properly.')
    sys.exit()
##/

## Move mash-tree ANI dendogram to main workdir and remove dRep-compare workdir
print('Retrieving PDF...')
cmd_cp = ['cp',workdir+'/glance_wd/figures/Primary_clustering_dendrogram.pdf',workdir+'/glance_clustering_dendogram.pdf']
subprocess.call(' '.join(map(str,cmd_cp)),shell=True)

if not keep_files:
    print('Cleaning workspace')
    shutil.rmtree(workdir+'/'+'glance_wd')
##/

## Print final
print('Output file: glance_clustering_dendogram.pdf')
##/