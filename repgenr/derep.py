#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import subprocess
from multiprocessing import Pool

### Script internal function
def round_down_to_nearest(number,round_to):
    return (number // round_to) * round_to
###/

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = 
'''
Cluster genomes based on average nucleotide identity (ANI) using the dRep software
''')

parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-sani','--secondary_ani',default=0.99,help='Average nucleotide identity (ANI) threshold for clustering in sensitive step of dRep (secondary clustering) (default: 0.99)')
parser.add_argument('-pani','--primary_ani',default=0.90,help='Average nucleotide identity (ANI) threshold for clustering in rough step of dRep (primary clustering) (default: 0.9)')
parser.add_argument('--pre_secondary_ani',default=None,help='On dataset chunks, ANI threshold for clustering in first stage of dRep. Only applied if -s/--process_size is >1 (default: same as -sani)')
parser.add_argument('--pre_primary_ani',default=None,help='On dataset chunks, ANI threshold for clustering in second stage of dRep. Only applied if -s/--process_size is >1 (default: same as -pani)')
parser.add_argument('--S_algorithm',default='fastANI',help='Algorithm to use in dRep fastANI. Possible values as of dRep v3.4.1: fastANI,gANI,goANI,ANIn,ANImf (default: fastANI)')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of total threads to use (default: 16)')
parser.add_argument('-p','--num_processes',type=int,default=1,help='Number of processes to execute in parallel\nMust be paired with -s/--process_size to split a large dataset into multiple batches that are processed separately.\In a second round, the dereplicated output of batches are subjected to dRep (default: 1)')
parser.add_argument('-s','--process_size',type=int,default=False,help='Number of genomes per process. Lower values increase speed (default: full dataset, i.e. no dataset splitting)')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
parser.add_argument('--skip_intermediary',action='store_true',help='If specified, will not attempt to recycle intermediary files from dRep')
parser.add_argument('--prompt',action='store_true',help='If specified and using the "chunking" feature (process_size < "number of input genomes") will prompt to continue to secondary dRep-run. This can be used to adjust parameters in the first round of dRep that is run on chunks.')
parser.add_argument('--secondary_only',action='store_true',help='If specified, will only do the second round of dRep on previously generated data from "chunks"')
parser.add_argument('--virus',action='store_true',help='If specified, will pass virus-tuned parameters to dRep.')
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

drep_secondary_run_only = args.secondary_only
drep_virus_mode = args.virus
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (make sure you have run "metadata" and "genome" commands):')
    print(workdir)
    sys.exit()
    
if drep_secondary_run_only and not os.path.exists(workdir+'/'+'derep_chunks_representative_genomes.tsv'):
    print('A secondary-only run was specified, but could not locate previous output files of stage 1! If you ran derep before, this may indicate that your dataset was processed as one batch in the first step - this makes the secondary step unnecessary')
    sys.exit()
#/
###/

### Check so genomes-dir exists (required by this script)
if not os.path.exists(workdir+'/'+'genomes'):
    print('Could not locate genome directory inside specified workdir (make sure you have run "metadata" and "genome" commands):')
    print(workdir)
    sys.exit()
###/

### WARNING: Put a warning for user if using pairwise aligners (e.g. ANImf as is run with --virus). It creates a squared number of intermediary files which leads to OS issues due to millions of files in a folder
num_genomes_files = len(os.listdir(workdir+'/'+'genomes'))
if (drep_virus_mode or S_algorithm != 'fastANI') and (num_genomes_files > 200 and (process_size == False or process_size > 200)):
    print('WARNING: a large number of genome-files are going to be run through dRep using a pairwise aligner. This generates a squared number of output that has multiple files per input, which may result in millions of files in a directory and subsequently slow down the OS.')
    print('For example, just the drep_workdir/log/cmd_logs/ will have 750\'000 files for 500 input sequences (3 files per sample, run in all-vs-all)')
    print('You can alleviate this issue by using --process_size <num_genomes_per_process> and it is recommended that this number is below 200')
    print('Recommended size for --process_size is between 100-200. A smaller number may lead to faster execution.')
    print('\nPausing for a few seconds. Did you read it yet? ;)',flush=True)
    # sleep for some seconds to user has a chance to see message
    import time
    time.sleep(5)
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
    print('A previous directory for drep chunks was found. It will be deleted now before creating the new one',flush=True)
    # For INFO, get number of files to remove. Found that sometimes clearing a previous workdir is EXTREMELY slow when the previous run generated alot of bullshit files
    print('Checking number of files to remove...',flush=True)
    num_files = 0
    warned_at = set()
    warned_at2 = set()
    for path,dirs,files in os.walk(chunks_workdir):
        num_files += len(files)
        round_500k = round_down_to_nearest(num_files,500000)
        if round_500k > 0 and not round_500k in warned_at:
            print('Current file count is '+'{:,} and counting...'.format(num_files),flush=True)
            warned_at.add(round_500k)
            
            round_1M = round_down_to_nearest(num_files,1000000)
            if round_1M > 0 and not round_1M in warned_at2:
                print('A large number of files can lead to long removal times. Please stay put',flush=True)
                warned_at2.add(round_1M)
    #/
    # Delete workding directory if it previously existed. We need to make sure the input to drep do not change if a user re-runs the command
    print('Number of files to remove: '+'{:,}'.format(num_files),flush=True)
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

#@ Check if there's a single genome in last chunk. If so, append it to next-to-last chunk
if downloaded_genomes_chunked and len(downloaded_genomes_chunked[-1]) == 1:
    print('When splitting input into chunks, the last chunk had only one genome. Will move this genome to another chunk.',flush=True)
    downloaded_genomes_chunked[-2] += downloaded_genomes_chunked[-1] # append last chunk to next-to-last
    del downloaded_genomes_chunked[-1] # remove last
#@/

if not drep_secondary_run_only: # dont print this status if we only do stage2
    print('Total of '+str(sum(map(len,downloaded_genomes_chunked)))+ ' downloaded genomes were split into '+str(len(downloaded_genomes_chunked))+' jobs',flush=True)
#/
# Soft-link genomes into chunked workspaces
print('Soft-linking genomes into subprocess workding directories',flush=True)
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
        print('Warning: Could not use intermediary files as they were empty or not located! Will need to re-calculate files',flush=True)
    else:
        print('Folder for intermediary files does not exist. Will create and populate this folder!',flush=True)
    skip_intermediary_readFiles = True
#/
#
if not skip_intermediary_readFiles:
    if os.path.exists(intermediary_files_dir) and os.path.exists(intermediary_files_dir+'/'+'Chdb.csv'):
        print('Found previously calculated files from dRep, will try to use those (dRep "Filtering" step)',flush=True)
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
                print('Warning: Did not find previous calculations for all genomes in chunk'+str(enum)+', will need to recompute these!',flush=True)
            #/
##/
###/

### Run drep on chunks (stage 1, dereplication intra-chunk)
## Define multiprocessing worker (used for stage1 and stage2)
def drep_worker(wd,pani,sani,threads,print_status):
    if print_status:            print(print_status + ' started...',flush=True)
    ## Run drep
    # compile base CMD
    drep_cmd = ['derep_worker.py','-wd',wd]
    #/
    # add virus-specific arguments to CMD
    if drep_virus_mode: # check if append virus-related settings (taken from: https://drep.readthedocs.io/en/latest/choosing_parameters.html#comparing-and-dereplicating-non-bacterial-genomes)
        drep_cmd += ['--S_algorithm','ANImf',
                     '-nc','0.5', # "cov_thresh"
                     '-N50W','0', # N50 weight
                     '-sizeW','1', # size weight
                     '--ignoreGenomeQuality',
                     '--clusterAlg','single']
    #/
    # add user arguments to CMD
    drep_cmd += ['--primary_ani',pani,'--secondary_ani',sani,'--threads',threads]
    if (not drep_virus_mode) or '--S_algorithm' in sys.argv: # add S_algorithm unless we are in virus-mode AND did not explicitly state S_algorithm
        drep_cmd += ['--S_algorithm',S_algorithm]
    #/
    # execute cmd
    subprocess.call(' '.join(map(str,drep_cmd)),shell=True)
    #/
    ##/
    # Check if output folder was created, else return false
    drep_output_dir = False
    if os.path.exists(wd+'/'+'genomes_derep_representants') and len(os.listdir(wd+'/'+'genomes_derep_representants')) != 0:
        drep_output_dir = wd
    #/
    if print_status:            print(print_status + ' finished!',flush=True)
    return drep_output_dir
##/
# Execute first round stage 1 (skip if user said to do round 2 only)
if not drep_secondary_run_only:
    ## setup multiprocessing
    # Check how many parallel processes can be run
    num_processes = min(num_processes,len(downloaded_genomes_chunked))
    num_threads_per_process = int(num_threads/num_processes)
    print('Processing genomes (number of parallel processes is '+str(num_processes)+ ' at '+str(num_threads_per_process)+' threads each)',flush=True)
    #/
    ##/
    ## start jobs
    pool = Pool(processes=num_processes)
    jobs = {}
    for enum,chunk_wd in enumerate(os.listdir(intra_chunks_wd)):
        chunk_wd_path = intra_chunks_wd+'/'+chunk_wd
        # Check so we have the expected folder structure
        if not os.path.exists(chunk_wd_path+'/'+'genomes'):
            print('Incorrect folder structure found for directory '+chunk_wd+'\n\tskipping...',flush=True)
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
    
    print('\tProcessing done!',flush=True)
    ##/
    ###/
    
    ### Fetch output of chunks (prepare for stage 2)
    print('Fetching genomes from stage 1',flush=True)
    # Save representative genomes and intermediary data from stage 1
    failed_jobs = []
    genomes_stage1 = []
    Chdb_genome_entries = {} # for intermediary file (Chdb)
    Chdb_header = 'header_placeholder\n' # will be re-assigned upon file reading
    nf = open(workdir+'/'+'derep_chunks_representative_genomes.tsv','w')
    for chunk,chunk_output_path in jobs_status.items():
        # Check if this chunk did not execute properly
        if chunk_output_path == False:
            print('Job '+chunk+' failed with error!',flush=True)
            failed_jobs.append(chunk)
            continue
        #/
        # write representative genomes
        for file_ in os.listdir(chunk_output_path+'/'+'genomes_derep_representants'):
            genomes_stage1.append(file_)
            
            writeArr = [file_,chunk]
            nf.write('\t'.join(map(str,writeArr))+'\n')
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
    nf.close()
    #/
    # Save Chdb intermediary file
    if not skip_intermediary:
        if not os.path.exists(intermediary_files_dir):      os.makedirs(intermediary_files_dir)
        with open(intermediary_files_dir+'/'+'Chdb.csv','w') as nf:
            nf.write(Chdb_header)
            for genome,row in Chdb_genome_entries.items():
                nf.write(row)
    #/
    # check if jobs failed
    if failed_jobs:
        print('Could not run dRep properly. Terminating.',flush=True)
        sys.exit()
    #/
    # print some info
    print('Number of genomes remaining after first round of dRep on chunks: '+str(len(genomes_stage1))+'/'+str(sum(map(len,downloaded_genomes_chunked))),flush=True)
    #/
#/
###/

### Run drep on chunks-output (stage 2, dereplication inter-chunk)
# Remove previous directory of genome representants
if os.path.exists(workdir+'/'+'genomes_derep_representants'):
    print('Previous directory of dereplicated genomes was erased',flush=True)
    shutil.rmtree(workdir+'/'+'genomes_derep_representants')
#/

drep_secondary_run = False
## Check if we had multiple chunks from stage 1 (then copy results to final folder)
if not drep_secondary_run_only and len(jobs_status) == 1:
    print('All datasets was processed in the same batch, will not do a second round of dRep',flush=True)
    chunk_output_path = list(jobs_status.items())[0][1]
    shutil.copytree(chunk_output_path+'/'+'genomes_derep_representants',inter_chunks_wd+'/'+'genomes_derep_representants')
    
    # remove the list of representative genomes (will re-create it from files in workdir/genomes_derep_representants)
    if os.path.exists(workdir+'/'+'derep_chunks_representative_genomes.tsv'):        os.remove(workdir+'/'+'derep_chunks_representative_genomes.tsv')
    #/
##/
## Else, run drep again
else:
    drep_secondary_run = True
    # Check if we prompt to continue before proceeding to 2nd round of dRep
    if prompt_before_step2:
        usr_inp = input('Proceed with second round of dRep? (no/n to terminate) ').strip()
        if usr_inp.lower() in ('n','no'):
            sys.exit('User input: '+usr_inp+', terminating software')
    #/
    # init genomes directory in inter-chunks (stage 2)
    inter_chunks_wd_genomes_path = inter_chunks_wd + '/' + 'genomes'
    if not os.path.exists(inter_chunks_wd_genomes_path):        os.makedirs(inter_chunks_wd_genomes_path)
    #/
    ## Setup folders and genomes
    # Pre: check if we run only stage2, then import some data from stage1
    if drep_secondary_run_only:
        genomes_stage1 = []
        with open(workdir+'/'+'derep_chunks_representative_genomes.tsv','r') as f:
            for line in f:
                line = line.strip('\n')
                genome,chunk = line.split('\t')
                genomes_stage1.append(genome)
    #/
    # Copy-in genomes from stage1 to stage2
    for genome in genomes_stage1:
        shutil.copy2(workdir+'/'+'genomes'+'/'+genome,inter_chunks_wd_genomes_path)
    #/
    ##/
    ## Setup intermediate data
    if not skip_intermediary:
        # init "drep_workdir" folder, and "data_tables" folder
        step2_drep_wd = inter_chunks_wd+'/'+'drep_workdir'
        step2_drep_data_tables = step2_drep_wd+'/'+'data_tables'
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
            print('Warning: Did not find previous calculations for all genomes to input for secondary dRep run (inter_chunks), will need to recompute these!',flush=True)
        #/
    ##/
    drep_worker(inter_chunks_wd,primary_ani,secondary_ani,num_threads,'[worker 1/1]')
##/
###/

### Fetch final dereplication output and tranfer to master-folder (main workingdirectory)
print('Finalizing dereplication',flush=True)
shutil.copytree(inter_chunks_wd+'/'+'genomes_derep_representants',workdir+'/'+'genomes_derep_representants')

# Write representative genomes
with open(workdir+'/'+'derep_representative_genomes.tsv','w') as nf:
    for file_ in os.listdir(workdir+'/'+'genomes_derep_representants'):
        nf.write(file_+'\n')
#/
###/
### Print number of genomes in dereplication-folder
print('Number of genomes remaining as representative genomes: '+str(len(os.listdir(workdir+'/'+'genomes_derep_representants')))+'/'+str(sum(map(len,downloaded_genomes_chunked))),flush=True)
###/

### Summarize compiled sequences
## Stage 1 (contained sequences during dRep of chunks)
if not drep_secondary_run_only:
    with open(workdir+'/'+'derep_chunks_clustered_genomes.tsv','w') as nf:
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

## Get genome stage1/stage2 statuses
genome_statuses_stage1 = {}
genome_statuses_stage2 = {}
# If we run secondary step only, then import old data
if drep_secondary_run_only:
    with open(workdir+'/'+'derep_genomes_status.tsv','r') as f:
        for line in f:
            line = line.strip('\n')
            genome,repgenome,stage1_status,stage2_status= line.split('\t')
            if stage1_status != 'NA':
                genome_statuses_stage1[genome] = stage1_status
            if stage2_status != 'NA':
                genome_statuses_stage1[genome] = stage1_status
#/
# 1. Get genome status from chunked workers (from intra_chunks wd, "stage 1"), skip this if we do stage 2 only
if not drep_secondary_run_only:
    for chunk,chunk_output_path in jobs_status.items():
        if chunk_output_path != False:
            with open(chunk_output_path+'/'+'derep_genome_status.tsv','r') as f:
                for line in f:
                    line = line.strip('\n')
                    genome,status = line.split('\t')
                    
                    genome_statuses_stage1[genome] = status
#/
# 2. Get genome status from inter_chunks wd "stage 2". Overwrite new information; since representative genomes from stage1 might be contained in stage2.
if drep_secondary_run:
    with open(inter_chunks_wd+'/'+'derep_genome_status.tsv','r') as f:
        for line in f:
            line = line.strip('\n')
            genome,status = line.split('\t')
            
            genome_statuses_stage2[genome] = status
#/
##/
# Mark which genomes were final representants
genomes_representants = list(os.listdir(workdir+'/'+'genomes_derep_representants'))
#/
## Write summary file (columns: file_name, representative (bool), QC status, ?
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
## Get genome quality data from dRep "genomeInformation.csv" files (genomes are removed by default in dRep if <75% ocmpleteness, >25% contamination, etc.)
# Parse data from work-folders
genomeInformation_header = None
genomeInformation = {} # genome name -> raw_row
for path,dirs,files in os.walk(chunks_workdir):
    for file_ in files:
        if file_ == 'genomeInformation.csv' and os.path.basename(path) == 'data_tables':
            with open(path+'/'+file_,'r') as f:
                for enum,line in enumerate(f):
                    line = line.strip('\n')
                    line = line.split(',')
                    # check if header-row, save it if we did not already
                    if enum == 0:
                        if genomeInformation_header == None:
                            genomeInformation_header = line
                        continue
                    #/
                    # Save rows
                    genome = line[0]
                    genomeInformation[genome] = line
                    #/
#/
# Check if there was a previous file, then parse data from it
if os.path.exists(workdir+'/'+'derep_genomeInformation.tsv'):
    with open(workdir+'/'+'derep_genomeInformation.tsv','r') as f:
        for enum,line in enumerate(f):
            line = line.strip('\n')
            line = line.split(',')
            # check if header-row, save it if we did not already
            if enum == 0:
                if genomeInformation_header == None:
                    genomeInformation_header = line
                continue
            #/
            # Save rows
            genome = line[0]
            genomeInformation[genome] = line
            #/
#/
# Dump file
with open(workdir+'/'+'derep_genomeInformation.tsv','w') as nf:
    rows_written = 0
    nf.write(','.join(genomeInformation_header)+'\n')
    for genome,row in genomeInformation.items():
        if 0 and not genome in genomes_all:
            print('Warning: genome had information stored, but was not located in "genomes"-folder. This indicates that this workdir has been re-run with new genomes. To be safe, previous derep-files should be removed if genome selection has changed.')
            continue
        nf.write(','.join(row)+'\n')
        rows_written += 1
#/
##/
## Output information of discarded files (using default genome filtering options described at "https://drep.readthedocs.io/en/latest/module_descriptions.html?highlight=GENOME%20FILTERING%20OPTION#dereplicate")
# dRep dereplicate:
# GENOME FILTERING OPTIONS:
#  -l LENGTH, --length LENGTH
#                        Minimum genome length (default: 50000)
#  -comp COMPLETENESS, --completeness COMPLETENESS
#                        Minumum genome completeness (default: 75)
#  -con CONTAMINATION, --contamination CONTAMINATION
#                        Maximum genome contamination (default: 25)
with open(workdir+'/'+'derep_potentially_discarded_genomes.tsv','w') as nf:
    for genome,row in genomeInformation.items():
        genome,completeness,contamination,strain_heterogeneity,length,N50,centrality = row
        try:        length = int(length)
        except:     length = 9999999999
        try:        completeness = float(completeness)
        except:     completeness = 999
        try:        contamination = float(contamination)
        except:     contamination = 999
        
        qc_fails = []
        if length < 5000:               qc_fails.append('fail_length')
        if completeness < 75:           qc_fails.append('fail_qc')
        if contamination > 25:          qc_fails.append('fail_contamination')
        
        if qc_fails:
            nf.write(genome+'\t'+','.join(qc_fails)+'\n')
##/
###/

### Clean up workspace
if not keep_files:
    print('Cleaning workspace',flush=True)
    # For INFO, get number of files to remove. Found that sometimes clearing a previous workdir is EXTREMELY slow when the previous run generated alot of bullshit files
    print('Checking number of files to remove...',flush=True)
    num_files = 0
    warned_at = set()
    warned_at2 = set()
    for path,dirs,files in os.walk(chunks_workdir):
        num_files += len(files)
        round_500k = round_down_to_nearest(num_files,500000)
        if round_500k > 0 and not round_500k in warned_at:
            print('Current file count is '+'{:,} and counting...'.format(num_files),flush=True)
            warned_at.add(round_500k)
            
            round_1M = round_down_to_nearest(num_files,1000000)
            if round_1M > 0 and not round_1M in warned_at2:
                print('A large number of files can lead to long removal times. Please stay put',flush=True)
                warned_at2.add(round_1M)
    
    print('Number of files to remove: '+'{:,}'.format(num_files),flush=True)
    #/
    shutil.rmtree(chunks_workdir)
###/

### Write parameters used
if not drep_secondary_run_only:
    with open(workdir+'/'+'derep_parameters.txt','w') as nf:
        nf.write('primary_ani'+'\t'+str(primary_ani)+'\n')
        nf.write('secondary_ani'+'\t'+str(secondary_ani)+'\n')
        nf.write('input_data_processed_in_two_stages'+'\t'+str(drep_secondary_run)+'\n')
        nf.write('number_of_dataset_splits'+'\t'+str(num_processes)+'\n')
        nf.write('datasets per split'+'\t'+str(process_size)+'\n')
        nf.write('pre_primary_ani'+'\t'+str(pre_primary_ani)+'\n')
        nf.write('pre_secondary_ani'+'\t'+str(pre_secondary_ani)+'\n')
        if (not drep_virus_mode) or '--S_algorithm' in sys.argv:        nf.write('S_algorithm'+'\t'+str(S_algorithm)+'\n') # write this only if not overridden by --virus or actually passed by the user
        nf.write('virus_mode'+'\t'+str(drep_virus_mode)+'\n')
else:
    # get old file
    old_lines = []
    with open(workdir+'/'+'derep_parameters.txt','r') as f:
        for line in f:
            if line.strip('\n') == '': break #check if empty line, means separator between stage1+stage2 run and stage2_only run
            old_lines.append(line)
    #/
    # write old file plus stage2_only
    with open(workdir+'/'+'derep_parameters.txt','w') as nf:
        nf.write(''.join(old_lines))
        
        nf.write('\n')
        nf.write('Secondary_only_parameters:\t\n')
        nf.write('primary_ani'+'\t'+str(primary_ani)+'\n')
        nf.write('secondary_ani'+'\t'+str(secondary_ani)+'\n')
        nf.write('S_algorithm'+'\t'+str(S_algorithm)+'\n')
    #/
###/

### Summarize contained/dereplicated-and-removed genomes.
# NOTE: I do this in a function so that previous uses of repgenr can summarize their derep results by [python -c "from derep import summarize_derep_genomes; workdir=<workdir path used when running repgenr>; summarize_derep_genomes()"]
# It can be moved into this script-file when previous use of repgenr has been barred
from derep_summarize import summarize_derep_genomes
summarize_derep_genomes(workdir)
###/