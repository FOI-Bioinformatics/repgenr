#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import subprocess
from multiprocessing import Pool

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-sani','--secondary_ani',default=0.99,help='Average nucleotide identity (ANI) threshold for clustering in sensitive step of dRep (secondary clustering)')
parser.add_argument('-pani','--primary_ani',default=0.90,help='Average nucleotide identity (ANI) threshold for clustering in rough step of dRep (primary clustering)')
parser.add_argument('--pre_secondary_ani',default=None,help='On dataset chunks, ANI threshold for clustering in first stage of dRep. Only applied if -s/--process_size is >1 (default: same as -sani)')
parser.add_argument('--pre_primary_ani',default=None,help='On dataset chunks, ANI threshold for clustering in second stage of dRep. Only applied if -s/--process_size is >1 (default: same as -pani)')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of threads to use (total)')
parser.add_argument('-p','--num_processes',type=int,default=1,help='Number of processes to execute in parallel\nMust be paired with -s/--process_size to split a large dataset into multiple batches that are processed separately. In a second round, the dereplicated output of batches are subjected to dRep (default: 1)')
parser.add_argument('-s','--process_size',type=int,default=False,help='Number of genomes per process. Lower values increase speed (default: full dataset, i.e. no dataset splitting)')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
#/
# parse input
args = parser.parse_args()
#args = parser.parse_args('-wd /data/users/jaclew/projects/2022/refDB_cond/tularensis_WIP -t 60 -p 1 -s 100'.split())

workdir = args.workdir

secondary_ani = args.secondary_ani
primary_ani = args.primary_ani
num_threads = args.threads
num_processes = args.num_processes
process_size = args.process_size

pre_secondary_ani = args.pre_secondary_ani
pre_primary_ani = args.pre_primary_ani
if pre_secondary_ani == None:       pre_secondary_ani = secondary_ani # set to same as secondary_ani if user did not specify value
if pre_primary_ani == None:         pre_primary_ani = primary_ani # set to same as primary_ani if user did not specify value

keep_files = args.keep_files
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (make sure you have run metadata and genome commands):')
    print(workdir)
    sys.exit()
#/
###/


### init workspace for chunking
# Define workingdir for chunking and make it (./derep_chunks_workdir with subfolder "intra_chunks" (stage1: drep on chunked genomes) and "inter_chunks" (stage2: drep on output from stage1)
# Tree structure:
# derep_chunk_workdir
# ├── intra_chunks (stage 1, optionally run as parallel processes)
# │   ├── chunk1
# │   │   ├── drep_out
# │   │   └── genomes
# │   ├── chunk2
# │   │   ├── drep_out
# │   │   └── genomes
# │   └── chunk3
# │       ├── drep_out
# │       └── genomes
# └── inter_chunks (stage 2)
#     ├── drep_out
#     └── genomes

chunks_workdir = workdir+'/'+'derep_chunks_workdir' # the results will be returned to derep_chunks_genomes
if os.path.exists(chunks_workdir):
    # Delete workding directory if it previously existed. We need to make sure the input to drep do not change if a user re-runs the command
    print('A previous directory for drep chunks was found. It will be deleted now before creating the new one')
    shutil.rmtree(chunks_workdir)
    #/
if not os.path.exists(chunks_workdir):     os.makedirs(chunks_workdir)

intra_chunks_wd = chunks_workdir+'/'+'intra_chunks'
inter_chunks_wd = chunks_workdir+'/'+'inter_chunks'
if not os.path.exists(intra_chunks_wd):     os.makedirs(intra_chunks_wd)
if not os.path.exists(inter_chunks_wd):     os.makedirs(inter_chunks_wd)
#/
# Divide downloaded genomes into chunks of user-specified size
downloaded_genomes_chunked = []
for file_ in os.listdir(workdir+'/'+'genomes'):
    if file_.endswith('.fasta') or file_.endswith('.fasta.gz'):
        # Check if init new chunk
        if not downloaded_genomes_chunked:                          downloaded_genomes_chunked.append([]) # no chunk was initiated
        if process_size != False and len(downloaded_genomes_chunked[-1]) > process_size:      downloaded_genomes_chunked.append([]) # previous chunk is filled
        #/
        downloaded_genomes_chunked[-1].append(file_)
        
print('Total of '+str(sum(map(len,downloaded_genomes_chunked)))+ ' downloaded genomes were split into '+str(len(downloaded_genomes_chunked))+' jobs')
#/
# Soft-link genomes into chunked workspaces
print('Soft-linking genomes into subprocess workding directories')
for enum,chunk in enumerate(downloaded_genomes_chunked):
    # init chunk wd
    chunk_wd = intra_chunks_wd + '/' + 'chunk'+str(enum)
    if not os.path.exists(chunk_wd):        os.makedirs(chunk_wd)
    #/
    # make genome-dir
    if not os.path.exists(chunk_wd+'/'+'genomes'):      os.makedirs(chunk_wd+'/'+'genomes')
    #/
    # softlink genomes
    for genome_fa in chunk:
        cmd_ln = ['ln','-s',os.path.abspath(workdir)+'/genomes/'+genome_fa,chunk_wd+'/'+'genomes']
        subprocess.call(' '.join(map(str,cmd_ln)),shell=True)
    #/
#/
###/

### Run drep on chunks (stage 1, dereplication intra-chunk)
## Define multiprocessing worker
def drep_worker(wd,pani,sani,threads):
    # Run drep
    drep_cmd = ['derep_worker.py','-wd',wd,'--primary_ani',pani,'--secondary_ani',sani,'--threads',threads]
    subprocess.call(' '.join(map(str,drep_cmd)),shell=True)
    #/
    # Check if output folder was created, else return false
    drep_output_dir = False
    if os.path.exists(chunk_wd+'/'+'genomes_derep_representants') and len(os.listdir(chunk_wd+'/'+'genomes_derep_representants')) != 0:
        drep_output_dir = chunk_wd
    #/
    return drep_output_dir
##/
## setup multiprocessing
# Check how many parallel processes can be run
num_processes = min(num_processes,len(downloaded_genomes_chunked))
num_threads_per_process = int(num_threads/num_processes)
print('Processing genomes (number of parallel processes is '+str(num_processes)+ ' at '+str(num_threads_per_process)+' threads each)')
##/
## start jobs
pool = Pool(processes=num_processes)
jobs = {}
for chunk_wd in os.listdir(intra_chunks_wd):
    chunk_wd_path = intra_chunks_wd+'/'+chunk_wd
    # Check so we have the expected folder structure
    if not os.path.exists(chunk_wd_path+'/'+'genomes'):
        print('Incorrect folder structure found for directory '+chunk_wd+'\n\tskipping...')
        continue
    #/
    jobs[chunk_wd] = pool.apply_async(drep_worker,args=(chunk_wd_path,pre_primary_ani,pre_secondary_ani,num_threads_per_process))
    #break
pool.close()
pool.join()
##/
## parse jobs
jobs_status = {}
for chunk_wd,job in jobs.items():
    jobs_status[chunk_wd] = job.get()

print('\tProcessing done!')
##/
###/

### Fetch output of chunks (prepare for stage 2)
print('Fetching genomes from stage 1')
# init genomes directory in inter-chunks (stage 2)
inter_chunks_wd_genomes_path = inter_chunks_wd + '/' + 'genomes'
if not os.path.exists(inter_chunks_wd_genomes_path):        os.makedirs(inter_chunks_wd_genomes_path)
#/
# Setup genome-directory in inter-chunks
genomes_stage1 = []
failed_jobs = []
for chunk,chunk_output_path in jobs_status.items():
    # Check if this chunk did not execute properly
    if chunk_output_path == False:
        print('Job '+chunk+' failed with error!')
        failed_jobs.append(chunk)
        continue
    #/
    
    # Fetch dereplicated genomes from chunk
    for file_ in os.listdir(chunk_output_path+'/'+'genomes_derep_representants'):
        genomes_stage1.append(file_)
        shutil.copy2(chunk_output_path+'/'+'genomes_derep_representants'+'/'+file_,inter_chunks_wd_genomes_path)
    #/
#/
if failed_jobs:
    print('Could not run dRep properly. Terminating.')
    sys.exit()
###/

### Run drep on chunks-output (stage 2, dereplication inter-chunk)
drep_secondary_run = False
## Check if we had multiple chunks from stage 1 (then copy results to final folder)
if len(jobs_status) == 1:
    print('All datasets was processed in the same batch, will not do a second round of dRep')
    chunk_output_path = list(jobs_status.items())[0][1]
    shutil.copytree(chunk_output_path+'/'+'genomes_derep_representants',inter_chunks_wd+'/'+'genomes_derep_representants',dirs_exist_ok=True)
##/
## Else, run drep again
else:
    print('Executing dRep with '+str(len(genomes_stage1))+ ' genomes')
    drep_secondary_run = True
    drep_worker(inter_chunks_wd,primary_ani,secondary_ani,num_threads)
##/
###/

### Fetch final dereplication output and tranfer to master-folder (main workingdirectory)
print('Finalizing dereplication')
shutil.copytree(inter_chunks_wd+'/'+'genomes_derep_representants',workdir+'/'+'genomes_derep_representants',dirs_exist_ok=True)

if not keep_files:
    print('Cleaning workspace')
    shutil.rmtree(chunks_workdir)
###/

### Get info about representative sequences and contained sequences
""" TO BE IMPLEMENTED. below-code is from derep_worker.py
## Parse genome representatives
genome_representatives = {} # representative_file_name -> cluster (cluster is parsed below)
for file_ in os.listdir(workdir+'/'+'drep_workdir/dereplicated_genomes'):
    genome_representatives[file_] = None
##/
## Write info about cluster-contained genomes
clusters_genomes = {}
with open(workdir+'/'+'drep_workdir/data_tables/Cdb.csv','r') as f:
    for enum,line in enumerate(f):
        if enum == 0: continue #skip header
        line = line.strip('\n')
        line = line.split(',')
        genome,clust = line[0],line[1]
        
        if not clust in clusters_genomes:        clusters_genomes[clust] = set()
        clusters_genomes[clust].add(genome)
        
        # check if current line is the representative
        if genome in genome_representatives:
            genome_representatives[genome] = clust
        #/
##/
## Write out
with open(workdir+'/'+'drep_clustered_genomes.tsv','w') as nf:
    for genome_rep,clust in genome_representatives.items():
        for genome_contained in clusters_genomes[clust]:
            if genome_contained == genome_rep: continue # skip self
            writeArr = [genome_rep,genome_contained]
            nf.write('\t'.join(writeArr)+'\n')
##/
###/
"""
### Write parameters used
with open(workdir+'/'+'derep_parameters.txt','w') as nf:
    nf.write('primary_ani'+'\t'+str(primary_ani)+'\n')
    nf.write('secondary_ani'+'\t'+str(secondary_ani)+'\n')
    nf.write('dRep was run twice'+'\t'+str(drep_secondary_run))
    nf.write('number of dataset splits'+'\t'+str(num_processes)+'\n')
    nf.write('datasets per split'+'\t'+str(process_size)+'\n')
    nf.write('pre_primary_ani'+'\t'+str(pre_primary_ani)+'\n')
    nf.write('pre_secondary_ani'+'\t'+str(pre_secondary_ani)+'\n')
###/