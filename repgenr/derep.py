#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import subprocess
from multiprocessing import Pool


"""
Example of chunking:
    Suppose we have 5000 genomes to dereplicate. Dereplicating them as one process takes too long time, therefore we wish to chunk the work up into 10 manageable batches.
    command: rengenr derep --process_size 500 --num_processes 10 --threads 50
    The command will divide the 5000 genomes into batches of 500 genomes each. It will run them all in parallel (as 10 different processes) at 50 threads total, 5 threads per process.
    If we want to make it less parallel, we lower the number of procsses
    command: rengenr derep --process_size 500 --num_processes 2 --threads 50
    Now it will run them 2 processses in parallel at 50 threads total, 25 threads per process.
"""
### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-sani','--secondary_ani',default=0.99,help='Average nucleotide identity (ANI) threshold for clustering in sensitive step of dRep (secondary clustering)')
parser.add_argument('-pani','--primary_ani',default=0.90,help='Average nucleotide identity (ANI) threshold for clustering in rough step of dRep (primary clustering)')
parser.add_argument('--pre_secondary_ani',default=None,help='On dataset chunks, ANI threshold for clustering in first stage of dRep. Only applied if -s/--process_size is >1 (default: same as -sani)')
parser.add_argument('--pre_primary_ani',default=None,help='On dataset chunks, ANI threshold for clustering in second stage of dRep. Only applied if -s/--process_size is >1 (default: same as -pani)')
parser.add_argument('--S_algorithm',default='fastANI',help='Algorithm to use in dRep fastANI. Possible values as of dRep v3.4.1: fastANI,gANI,goANI,ANIn,ANImf (default: fastANI)')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of total threads to use (default: 16)')
parser.add_argument('-p','--num_processes',type=int,default=1,help='Number of processes to execute in parallel\nMust be paired with -s/--process_size to split a large dataset into multiple batches that are processed separately.\In a second round, the dereplicated output of batches are subjected to dRep (default: 1)')
parser.add_argument('-s','--process_size',type=int,default=False,help='Number of genomes per process. Lower values increase speed (default: full dataset, i.e. no dataset splitting)')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
parser.add_argument('--skip_intermediary',action='store_true',help='If specified, will not attempt to recycle intermediary files from dRep')
parser.add_argument('--prompt',action='store_true',help='If specified and using the "chunking" feature (process_size < "number of input genomes") will prompt to continue to secondary dRep-run. This can be used to adjust parameters in the first round of dRep that is run on chunks.')
#/
# parse input
args = parser.parse_args()
#args = parser.parse_args('-wd /data/users/jaclew/projects/2022/refDB_cond/tularensis_WIP_sub50 -t 20 -p 1 -s 100 --keep_files'.split())

workdir = args.workdir

secondary_ani = args.secondary_ani
primary_ani = args.primary_ani
S_algorithm = args.S_algorithm
num_threads = args.threads
num_processes = args.num_processes
process_size = args.process_size

pre_secondary_ani = args.pre_secondary_ani
pre_primary_ani = args.pre_primary_ani
if pre_secondary_ani == None:       pre_secondary_ani = secondary_ani # set to same as secondary_ani if user did not specify value
if pre_primary_ani == None:         pre_primary_ani = primary_ani # set to same as primary_ani if user did not specify value

keep_files = args.keep_files
skip_intermediary = args.skip_intermediary
prompt_before_step2 = args.prompt
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (make sure you have run "metadata" and "genome" commands):')
    print(workdir)
    sys.exit()
#/
###/


### init workspace for chunking
# Define workingdir for chunking and make it (./dereplication_workdir with subfolder "intra_chunks" (stage1: drep on chunked genomes) and "inter_chunks" (stage2: drep on output from stage1)
# Tree structure:
# dereplication_workdir
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

chunks_workdir = workdir+'/'+'dereplication_workdir' # the results will be returned to derep_chunks_genomes
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
        if process_size != False and len(downloaded_genomes_chunked[-1]) >= process_size:      downloaded_genomes_chunked.append([]) # previous chunk is filled
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
## Check if we want to use intermediary dRep data (Chdb.csv)
# Check so that the intermediary file exists
intermediary_files_dir = workdir+'/'+'dereplication_intermediary_files'
skip_intermediary_readFiles = False
if (not skip_intermediary) and (not os.path.exists(intermediary_files_dir+'/'+'Chdb.csv') or os.path.getsize(intermediary_files_dir+'/'+'Chdb.csv') == 0):
    if os.path.exists(intermediary_files_dir):
        print('Warning: Could not use intermediary files as they were empty or not located! Will need to re-calculate files')
    else:
        print('Folder for intermediary files does not exist. Will create and populate this folder!')
    skip_intermediary_readFiles = True
#/
#
if not skip_intermediary_readFiles:
    if os.path.exists(intermediary_files_dir) and os.path.exists(intermediary_files_dir+'/'+'Chdb.csv'):
        print('Found previously calculated files from dRep, will try to use those (dRep "Filtering" step)')
        for enum,chunk in enumerate(downloaded_genomes_chunked):
            chunk_wd_path = intra_chunks_wd + '/' + 'chunk'+str(enum)
            
            # init "drep_workdir" folder and "data_tables" folder
            chunk_drep_wd = chunk_wd_path+'/'+'drep_workdir'
            chunk_drep_data_tables = chunk_wd_path+'/'+'drep_workdir'+'/'+'data_tables'
            if not os.path.exists(chunk_drep_wd):                        os.makedirs(chunk_drep_wd)
            if not os.path.exists(chunk_drep_data_tables):      os.makedirs(chunk_drep_data_tables)
            #/
            # Get genomes of chunk from intermediary files
            chunk_genome_entries = {} # genome -> Chdb row
            Chdb_header = 'header_placeholder\n' # set default
            with open(intermediary_files_dir+'/'+'Chdb.csv','r') as f:
                for ln,line_raw in enumerate(f):
                    # parse header
                    if ln == 0:
                        Chdb_header = line_raw
                        continue
                    #/
                    # parse row
                    line = line_raw # keep raw line for instant write
                    line = line.strip('\n')
                    line = line.split(',')
                    
                    genome = line[0]
                    if genome in chunk:
                        chunk_genome_entries[genome] = line_raw
                    #/
            #/
            # Check if we found all genomes for chunk in Chdb, then we can write this file
            if set(chunk_genome_entries) == set(chunk):
                with open(chunk_drep_data_tables+'/'+'Chdb.csv','w') as nf:
                    # write header
                    nf.write(Chdb_header)
                    #/
                    # write rows
                    for genome,entry in chunk_genome_entries.items():
                        nf.write(entry)
                    #/
            else:
                print('Warning: Did not find previous calculations for all genomes in chunk'+str(enum)+', will need to recompute these!')
            #/
##/
###/

### Run drep on chunks (stage 1, dereplication intra-chunk)
## Define multiprocessing worker
def drep_worker(wd,pani,sani,threads,print_status):
    if print_status:            print(print_status + ' started...')
    # Run drep
    drep_cmd = ['derep_worker.py','-wd',wd,'--primary_ani',pani,'--secondary_ani',sani,'--threads',threads,'--S_algorithm',S_algorithm]
    subprocess.call(' '.join(map(str,drep_cmd)),shell=True)
    #/
    # Check if output folder was created, else return false
    drep_output_dir = False
    if os.path.exists(wd+'/'+'genomes_derep_representants') and len(os.listdir(wd+'/'+'genomes_derep_representants')) != 0:
        drep_output_dir = wd
    #/
    if print_status:            print(print_status + ' finished!')
    return drep_output_dir
##/
## setup multiprocessing
# Check how many parallel processes can be run
num_processes = min(num_processes,len(downloaded_genomes_chunked))
num_threads_per_process = int(num_threads/num_processes)
print('Processing genomes (number of parallel processes is '+str(num_processes)+ ' at '+str(num_threads_per_process)+' threads each)')
#/
##/
## start jobs
pool = Pool(processes=num_processes)
jobs = {}
for enum,chunk_wd in enumerate(os.listdir(intra_chunks_wd)):
    chunk_wd_path = intra_chunks_wd+'/'+chunk_wd
    # Check so we have the expected folder structure
    if not os.path.exists(chunk_wd_path+'/'+'genomes'):
        print('Incorrect folder structure found for directory '+chunk_wd+'\n\tskipping...')
        continue
    #/
    print_status = '[worker '+str(enum+1)+'/'+str(len(downloaded_genomes_chunked))+']' # Idea: Print-status holds the ID of the worker. It will display " [worker 1/50] " (input ID) followed by "started..." and --::-- " finished!" (appended in worker-function)
    jobs[chunk_wd] = pool.apply_async(drep_worker,args=(chunk_wd_path,pre_primary_ani,pre_secondary_ani,num_threads_per_process,print_status))
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
Chdb_genome_entries = {} # for intermediary file (Chdb)
Chdb_header = 'header_placeholder\n' # will be re-assigned upon file reading
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
    # Fetch intermediary data that we will save to outer
    if not skip_intermediary:
        if os.path.exists(chunk_output_path+'/'+'drep_workdir'+'/'+'data_tables'+'/'+'Chdb.csv'):
            with open(chunk_output_path+'/'+'drep_workdir'+'/'+'data_tables'+'/'+'Chdb.csv','r') as f:
                for ln,line_raw in enumerate(f):
                    # parse header
                    if ln == 0:
                        Chdb_header = line_raw
                        continue
                    #/
                    # parse row
                    line = line_raw # keep raw line for instant write
                    line = line.strip('\n')
                    line = line.split(',')
                    
                    genome = line[0]
                    Chdb_genome_entries[genome] = line_raw
                    #/
    #/
#/
# Save Chdb intermediary file
if not skip_intermediary:
    if not os.path.exists(intermediary_files_dir):      os.makedirs(intermediary_files_dir)
    with open(intermediary_files_dir+'/'+'Chdb.csv','w') as nf:
        nf.write(Chdb_header)
        for genome,row in Chdb_genome_entries.items():
            nf.write(row)
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
    print('Number of genomes remaining after first round of dRep on chunks: '+str(len(genomes_stage1))+'/'+str(sum(map(len,downloaded_genomes_chunked))))
    drep_secondary_run = True
    # Check if we prompt to continue before proceeding to 2nd round of dRep
    if prompt_before_step2:
        usr_inp = input('Continue with second round of dRep? (y/n) ').strip()
        if not usr_inp.lower() in ('y','yes'):
            sys.exit('User input: '+usr_inp+', terminating sotware')
    #/
    ## Setup intermediate data
    if not skip_intermediary:
        # init "drep_workdir" folder and "data_tables" folder
        step2_drep_wd = inter_chunks_wd+'/'+'drep_workdir'
        step2_drep_data_tables = inter_chunks_wd+'/'+'drep_workdir'+'/'+'data_tables'
        if not os.path.exists(step2_drep_wd):                        os.makedirs(step2_drep_wd)
        if not os.path.exists(step2_drep_data_tables):      os.makedirs(step2_drep_data_tables)
        #/
        # Get genomes of step1 (dereplicated genomes in step1) from intermediary files
        stage1_genome_entries = {} # genome -> Chdb row
        Chdb_header = 'header_placeholder\n' # set default
        with open(intermediary_files_dir+'/'+'Chdb.csv','r') as f:
            for ln,line_raw in enumerate(f):
                # parse header
                if ln == 0:
                    Chdb_header = line_raw
                    continue
                #/
                # parse row
                line = line_raw # keep raw line for instant write
                line = line.strip('\n')
                line = line.split(',')
                
                genome = line[0]
                if genome in genomes_stage1:
                    stage1_genome_entries[genome] = line_raw
                #/
        #/
        # Check if we found all genomes for input to stage2 in Chdb, then we can write this file
        if set(stage1_genome_entries) == set(genomes_stage1):
            with open(step2_drep_data_tables+'/'+'Chdb.csv','w') as nf:
                # write header
                nf.write(Chdb_header)
                #/
                # write rows
                for genome,entry in stage1_genome_entries.items():
                    nf.write(entry)
                #/
        else:
            print('Warning: Did not find previous calculations for all genomes to input for secondary dRep run (inter_chunks), will need to recompute these!')
        #/
    ##/
    drep_worker(inter_chunks_wd,primary_ani,secondary_ani,num_threads,'[worker 1/1]')
##/
###/

### Fetch final dereplication output and tranfer to master-folder (main workingdirectory)
print('Finalizing dereplication')
shutil.copytree(inter_chunks_wd+'/'+'genomes_derep_representants',workdir+'/'+'genomes_derep_representants',dirs_exist_ok=True)
###/
### Print number of genomes in dereplication-folder
print('Number of genomes remaining as representative genomes: '+str(len(os.listdir(workdir+'/'+'genomes_derep_representants')))+'/'+str(sum(map(len,downloaded_genomes_chunked))))
###/

### Summarize compiled sequences
## Stage 1 (contained sequences during dRep of chunks)
with open(workdir+'/'+'derep_chunks_clustered_genomes.tsv','w') as nf:
    if 0 and 'write header?':
        writeArr = ['representant','dereplicated'] # write header
        if drep_secondary_run:  writeArr.append('chunk')    # append column if we ran dRep in two stages. Stage 1 shall include a column saying which chunk it was processed in.
        nf.write('\t'.join(writeArr)+'\n')
        
    for chunk,chunk_output_path in jobs_status.items():
        if chunk_output_path != False:
            # Get genome representative<->contained relation
            with open(chunk_output_path+'/'+'derep_clustered_genomes.tsv','r') as f:
                for line in f:
                    line = line.strip('\n')
                    writeArr = line.split('\t')
                    if drep_secondary_run:  writeArr.append(chunk) # if we ran dRep in two stages, write representant + chunk from file, and add info about which chunk this relation was processed from
                    nf.write('\t'.join(writeArr)+'\n')
            #/
##/
## Stage 2 (contained sequences during dRep of stage1 output)
# Check if stage2 was not run (all datasets were run as one process in stage1)
if not drep_secondary_run:
    # Rename file from stage 1
    os.rename(workdir+'/'+'derep_chunks_clustered_genomes.tsv',workdir+'/'+'derep_clustered_genomes.tsv')
    #/
#/
# Else, copy file from stage 2
else:
    shutil.copy2(inter_chunks_wd+'/'+'derep_clustered_genomes.tsv',workdir+'/'+'derep_clustered_genomes.tsv')
#/
##/
###/

### Summarize dataset statuses
## Get genomes
genomes_all = list(os.listdir(workdir+'/'+'genomes'))
##/
## 1. Get genome status from chunked workers (from intra_chunks wd, "stage 1")
genome_statuses_stage1 = {}
for chunk,chunk_output_path in jobs_status.items():
    if chunk_output_path != False:
        with open(chunk_output_path+'/'+'derep_genome_status.tsv','r') as f:
            for line in f:
                line = line.strip('\n')
                genome,status = line.split('\t')
                
                genome_statuses_stage1[genome] = status
##/
## 2. Get genome status from inter_chunks wd "stage 2". Overwrite new information; since representative genomes from stage1 might be contained in stage2.
genome_statuses_stage2 = {}
if drep_secondary_run:
    with open(inter_chunks_wd+'/'+'derep_genome_status.tsv','r') as f:
        for line in f:
            line = line.strip('\n')
            genome,status = line.split('\t')
            
            genome_statuses_stage2[genome] = status
##/
## Mark which genomes were final representants
genomes_representants = list(os.listdir(workdir+'/'+'genomes_derep_representants'))
##/
## Write summary file
with open(workdir+'/'+'derep_genomes_status.tsv','w') as nf:
    # Header:
    # genome, representative genome (bool), stage1 status, stage2 status
    for genome in sorted(genomes_all):
        writeArr = []
        # genome
        writeArr.append(genome)
        #/
        # is_representative?
        writeArr.append(genome in genomes_representants) # true/false
        #/
        # stage1 status
        if genome in genome_statuses_stage1:
            writeArr.append(genome_statuses_stage1[genome])
        else:
            writeArr.append('NA')
        #/
        # stage2 status
        if genome in genome_statuses_stage2:
            writeArr.append(genome_statuses_stage2[genome])
        else:
            writeArr.append('NA')
        #/
        # Write
        nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
##/
###/

### Clean up workspace
if not keep_files:
    print('Cleaning workspace')
    shutil.rmtree(chunks_workdir)
###/

### Write parameters used
with open(workdir+'/'+'derep_parameters.txt','w') as nf:
    nf.write('primary_ani'+'\t'+str(primary_ani)+'\n')
    nf.write('secondary_ani'+'\t'+str(secondary_ani)+'\n')
    nf.write('input_data_processed_in_two_stages'+'\t'+str(drep_secondary_run)+'\n')
    nf.write('number_of_dataset_splits'+'\t'+str(num_processes)+'\n')
    nf.write('datasets per split'+'\t'+str(process_size)+'\n')
    nf.write('pre_primary_ani'+'\t'+str(pre_primary_ani)+'\n')
    nf.write('pre_secondary_ani'+'\t'+str(pre_secondary_ani)+'\n')
    nf.write('S_algorithm'+'\t'+str(S_algorithm)+'\n')
###/