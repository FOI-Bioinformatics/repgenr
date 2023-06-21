#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess
import argparse
import gzip
import ast
from multiprocessing import Pool

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m','--mode',required=True,choices=('accurate','fast',),help='Strategy to use for tree generation (accurate: pairwise progressivemauve+iqtree, fast: mashtree)')
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of total threads to use (default: 16)')
parser.add_argument('--no_outgroup',action='store_true',help='If specified, will not include the outgroup organism into the tree calculation')
parser.add_argument('--all_genomes',action='store_true',help='If specified, will run on all genomes and not on de-replicated genomes')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files (accruate-mode only)')
#/
# parse input
args = parser.parse_args()
num_threads = args.threads
run_mode = args.mode

workdir = args.workdir
skip_outgroup = args.no_outgroup
run_on_all_genomes = args.all_genomes

keep_files = args.keep_files
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

#### Generate tree using progressivemauve+x2fa/iqtree
if run_mode == 'accurate':
    ### Init phylo workdir (progressivemauve + x2fa and iqtree produce output files)
    phylo_wd = workdir+'/'+'phylo_workdir'
    # Remove previous phylo wd if it exists
    if os.path.exists(phylo_wd):
        shutil.rmtree(phylo_wd)
    #/
    # init phy wd
    os.makedirs(phylo_wd)
    #/
    ###/
    
    ### 
    ## Run Progressivemauve (pairwise)
    # find genomes -name "*.fasta" -exec basename {} \; | xargs -P 40 -I {} sh -c 'progressiveMauve --output progressivemauve/{}.xmfa outgroup/Rhizobiaceae_Ochrobactrum_B_teleogrylli_GCF_006376685.1.fasta genomes/{}'
    ## Run x2fa (xmfa to fa)
    # find progressivemauve/ -name "*.xmfa" -exec basename {} \; | xargs -I {} sh -c '../x2fa.py progressivemauve/{} outgroup/Rhizobiaceae_Ochrobactrum_B_teleogrylli_GCF_006376685.1.fasta 0 progressivemauve/{}.fa'
    ## Merge fa to multi-fa
    # parse .fa files and extract alignment for "rname" (first file only) and then "qname"
    
    ### Run progressivemauve and convert to fasta
    pmauve_ref = None
    pmauve_genome_list = []
    ## Find genomes to use. Softlink to new directory (progressivemauve produes additional files at input location...)
    # Get genomes
    for enum,file_ in enumerate(sorted(os.listdir(genome_files_dir))):
        # make first file the reference
        if enum == 0:
            pmauve_ref = file_
        #/
        pmauve_genome_list.append(genome_files_dir+'/'+file_)
    #/
    # Check if add in outgroup
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
        # Add outgroup to genome-list
        pmauve_genome_list.append(workdir+'/'+'outgroup'+'/'+outgroup_file)
        #/
    #/
    # soft-link genomes
    print('Soft-linking genomes...')
    pmauve_genomes_dir = phylo_wd+'/'+'pmauve_genomes'
    os.makedirs(pmauve_genomes_dir)
    
    for genome_file in pmauve_genome_list:
        genome_basename = os.path.basename(genome_file)
        cmd_ln = ['ln','-s',genome_file,pmauve_genomes_dir+'/'+genome_basename]
        subprocess.call(' '.join(map(str,cmd_ln)),shell=True)
    #/
    ##/
    ## Execute progressivemauve and convert to fasta
    # Validate that we have a reference
    if not pmauve_ref:
        print('Was unable to determine a reference to use in alignment! Make sure you have genomes...')
        sys.exit()
    #/
    
    # Define progressivemauve+x2fa worker
    def pmauve_x2fa_worker(output_path,reference_path,query_path,print_status=''):
        if print_status:            print(print_status + ' started...')
        
        query_name = os.path.splitext(query_path.split('/')[-1])[0]
        query_xmfa_out = output_path+'/'+query_name+'.xmfa'
        query_fa_out = output_path+'/'+query_name+'.fa'
        
        # Run progressiveMauve
        print('\nExecuting progressiveMauve')
        pmauve_cmd = ['progressiveMauve','--output',query_xmfa_out,reference_path,query_path]
        subprocess.call(' '.join(map(str,pmauve_cmd)),shell=True)
        #/
        
        # Run x2fa
        print('\nExecuting x2fa.py')
        pmauve_cmd = ['/data/users/jaclew/projects/2022/refDB_cond/x2fa.py',query_xmfa_out,reference_path,0,query_fa_out]
        subprocess.call(' '.join(map(str,pmauve_cmd)),shell=True)
        #/
        
        # Check if output file was created, else return false
        have_output = False
        if os.path.exists(query_fa_out) and os.path.getsize(query_fa_out) > 0:
            pmauve_has_output = True
        #/
        if print_status:            print(print_status + ' finished!')
        return pmauve_has_output
    #/
    # setup multiprocessing
    pmauve_outdir = phylo_wd+'/'+'progressivemauve'
    os.makedirs(pmauve_outdir)
    print('Processing genomes (number of parallel processes is '+str(num_threads)+ ' at '+str(1)+' threads each)') # progressiveMauve runs single-threaded
    #/
    # start jobs
    pool = Pool(processes=num_threads)
    jobs = {}
    for enum,genome_file in enumerate(os.listdir(pmauve_genomes_dir)):
        if genome_file == pmauve_ref: continue # skip self-alignment of reference
        print_status = '[worker '+str(enum+1)+'/'+str(len(pmauve_genome_list))+']' # Idea: Print-status holds the ID of the worker. It will display " [worker 1/50] " (input ID) followed by "started..." and --::-- " finished!" (appended in worker-function)
        
        jobs[enum] = pool.apply_async(pmauve_x2fa_worker,args=(pmauve_outdir,pmauve_genomes_dir+'/'+pmauve_ref,pmauve_genomes_dir+'/'+genome_file,print_status))
        #break
    pool.close()
    pool.join()
    #/
    # parse jobs
    jobs_status = {}
    for chunk_wd,job in jobs.items():
        jobs_status[chunk_wd] = job.get()

    print('\tProcessing done!')
    
    if not all(list(jobs_status.values())):
        print('Some jobs failed. Terminating.')
        sys.exit()
    #/
    ##/
    ### Concatenate single-alignment-fasta to multi-alignment-fasta and fix fasta-headers
    print('Concatenating alignment-fasta files to MSA-fasta...')
    ref_written = False
    write_line = True
    with open(phylo_wd+'/'+'msa.fasta','w') as nf:
        for file_ in os.listdir(pmauve_outdir):
            if file_.endswith('.fa'):
                with open(pmauve_outdir+'/'+file_,'r') as f:
                    for line in f:
                        # Check if header-line. Then remove path to file and keep only the file (genome) name.
                        if line[0] == '>':
                            line = '>' + os.path.splitext(line.split('/')[-1])[0]+'\n' # keeps only the file name
                        #/
                        
                        # Check if toggle writing (the "ref" exist in all single-alignment-fastas, but we just want one entry for the MSA-fasta)
                        if line[0] == '>':
                            if line.find(os.path.splitext(pmauve_ref)[0]) != -1:
                                if ref_written: # if ref was written, dont write more ref-lines
                                    write_line = False
                                else: # Write ref once
                                    write_line = True
                                    ref_written = True
                            else:
                                write_line = True
                        #/
                        
                        # write
                        if write_line:
                            nf.write(line)
                        #/
    ###/
    
    ### Run IQTREE (inside phylo_wd, iqtree produces some file)
    cmd_iqtree = ['iqtree','-T','auto','--threads-max',num_threads,'-s',phylo_wd+'/'+'msa.fasta']
    # check if using outgroup
    if not skip_outgroup:
        cmd_iqtree += ['-o',outgroup_file.replace('.fasta','')]
    #/
    subprocess.call(' '.join(map(str,cmd_iqtree)),shell=True)
    ###/
    
    ### Save which sample was used as reference in Progressivemauve and copy-out tree-file
    with open(phylo_wd+'/'+'tree_reference.txt','w') as nf:
        nf.write(pmauve_ref+'\n')
    
    with open(phylo_wd+'/'+'msa.fasta.treefile','r') as f:
        with open(output_file_path,'w') as nf:
            for line in f:
                nf.write(line)
    ###/
    ### Clean up workspace
    if not keep_files:
        print('Cleaning workspace')
        shutil.rmtree(phylo_wd)
    ###/
####/

#### Generate tree using mashtree
if run_mode == 'fast':
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
####/
