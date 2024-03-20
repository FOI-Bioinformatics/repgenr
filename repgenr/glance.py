#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
from matplotlib import pyplot as plt


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = 
'''
Perform rapid overview of ANI-clusters
''')

parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-t','--threads',type=int,default=24,help='Number of threads to use')
parser.add_argument('--plot_max',type=float,default=1,help='Maximum value to be included in boxplot and histogram (default:not set)')
parser.add_argument('--plot_min',type=float,default=0,help='Minimum value to be included in boxplot and histogram (default:not set)')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
#/
# parse input
args = parser.parse_args()

workdir = args.workdir
num_threads = args.threads

plot_max_val = args.plot_max
plot_min_val = args.plot_min

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

### Check if previous workdir of glace exist. Then remove it, since dRep otherwise complains about "not not provide genomes"
if os.path.exists(workdir+'/glance_wd'):
    shutil.rmtree(workdir+'/glance_wd')
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
##/

## Compute a boxplot of ANI values from the Mdb.csv file
print('Computing boxplot and histogram...')
if os.path.exists(workdir+'/glance_wd/data_tables/Mdb.csv'):
    # Parse similarity values from file
    similarity_values = []
    with open(workdir+'/glance_wd/data_tables/Mdb.csv','r') as f:
        for ln,line in enumerate(f):
            if ln == 0: continue # skip header line
            line = line.strip('\n')
            line = line.split(',')
            genome1,genome2,dist,similarity = line
            try:
                similarity = float(similarity)
            except:
                print('[WARNING] could not parse similarity-value as a float: '+similarity+', skipping...')
            if genome1 == genome2: continue # skip if current comparison is to self
            
            # check if user specified cutoffs or thresholds
            if similarity > plot_max_val: continue
            if similarity < plot_min_val: continue
            #/
            
            similarity_values.append(similarity)
    #/
    
    # Make plot
    fig, ax = plt.subplots()
    ax.boxplot(similarity_values)
    ax.set_ylabel('MASH ANI')
    ax.set_title('MASH ANI all-vs-all comparison. Number of values: '+str(len(similarity_values)))
    plt.tight_layout()
    plt.savefig(workdir+'/glance_MASH_ANI_similarity_boxplot.png')
    
    fig, ax = plt.subplots()
    ax.hist(similarity_values,bins=100)
    ax.set_ylabel('MASH ANI')
    ax.set_title('MASH ANI all-vs-all comparison. Number of values: '+str(len(similarity_values)))
    plt.tight_layout()
    plt.savefig(workdir+'/glance_MASH_ANI_similarity_histogram.png')
    #/
else:
    print('Could not locate Mdb.csv file. Boxplot/Histogram files will not be output. Make sure dRep compare was run properly')
##/

## Clean workspace
if not keep_files:
    print('Cleaning workspace')
    shutil.rmtree(workdir+'/'+'glance_wd')
##/

## Print final
print('Output files: glance_clustering_dendogram.pdf, glance_MASH_ANI_similarity_boxplot.png, glance_MASH_ANI_similarity_histogram.png')
##/