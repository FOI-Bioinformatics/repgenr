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
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = 
'''
Compute phylogenetic tree based on all or representative genomes
''')


parser.add_argument('-m','--mode',required=True,choices=('accurate','fast',),help='Strategy to use for tree generation (accurate: pairwise progressivemauve+iqtree, fast: mashtree)')
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of total threads to use (default: 16)')
parser.add_argument('-B','--bootstrap',type=int,default=0,help='Number of bootstrap iterations. Only applicable to "--mode accurate". Passes value to IQ-TREE "ultrafast bootstrap" (iqtree -B parameter). Minimum value 1000. (default: 0)')
parser.add_argument('--no_outgroup',action='store_true',help='If specified, will not include the outgroup organism into the tree calculation')
parser.add_argument('--all_genomes',action='store_true',help='If specified, will run on all genomes and not on de-replicated genomes')
parser.add_argument('--halt_after_msa',action='store_true',help='If specified, will stop execution the multiple sequence alignment file is generated. Implies "--keep_msa". This argument is useful when RepGenR is used as part of another pipe (accurate-mode only)')
parser.add_argument('--keep_msa',action='store_true',help='If specified, will save multiple sequence alignment file "msa_[derep/all].fasta" (accurate-mode only)')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files (accurate-mode only)')
parser.add_argument('--progressivemauve_ref',required=False,default=None,help='Accurate mode only. Specify file-name including extension of reference to use during progressivemauve (default: use first genome file in directory)')
parser.add_argument('--progressivemauve_single',action='store_true',help='If specified, will not parallelize progressivemauve. Sometimes this leads to a yet unknown error. Running one process at a time has solved the issue')
#/
# parse input
args = parser.parse_args()
num_threads = args.threads
run_mode = args.mode

bootstrap_num_iterations = args.bootstrap

workdir = args.workdir
skip_outgroup = args.no_outgroup
run_on_all_genomes = args.all_genomes

halt_after_msa = args.halt_after_msa

keep_msa = args.keep_msa
keep_files = args.keep_files

pmauve_ref = args.progressivemauve_ref
progressivemauve_single = args.progressivemauve_single
#/
## validate input
# check dirs
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
# check bootstrap input
if bootstrap_num_iterations > 0 and bootstrap_num_iterations < 1000:
    print('Bootstrap value must be at least 1000. For more details, see IQ-TREE documentation, section \n\tULTRAFAST BOOTSTRAP/JACKKNIFE: -B \t Replicates for ultrafast bootstrap')
    sys.exit()
#/
# Check keep_msa/halt_after_msa (requires --mode accurate)
if (keep_msa or halt_after_msa) and not run_mode == 'accurate':
    print('Arugment --keep_msa (including --halt_after_msa) requires --run_mode=accurate to run\nTerminating!')
    sys.exit()
#/
# check halt_after_msa (implies --keep_msa; toggle on --keep_msa if it was not done)
if halt_after_msa and not keep_msa:
    print('[INFO] --keep_msa was not specified: Setting --keep_msa to True')
    keep_msa = True
#/
##/
###/

### Get absolute path for workdir, needed for software calls
exec_cwd = os.getcwd()
if not workdir[0] == '/':       workdir = exec_cwd + '/' + workdir # get absolute path to workdir unless it was already supplied (full paths start with a slash, i.e. "/path/to/my/wd" as opposed to relative path "wd/")
###/

### Set input/output paths depending on input type (all genomes vs. dereplicated genomes)
# Check if run on de-rep genomes
if not run_on_all_genomes:
    genome_files_dir = workdir+'/'+'genomes_derep_representants'
    output_file_path = workdir+'/'+'genomes_derep_representants.dnd'
    msa_output_file_path = workdir+'/'+'msa_derep.fasta' # for --mode accurate only
#/
# Else, run on all genomes
else:
    genome_files_dir = workdir+'/'+'genomes'
    output_file_path = workdir+'/'+'genomes.dnd'
    msa_output_file_path = workdir+'/'+'msa_all.fasta' # for --mode accurate only
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
    
    ### Run progressivemauve and convert to fasta
    pmauve_genome_list = []
    ## Find genomes to use. Softlink to new directory (progressivemauve produes additional files at input location...)
    # Get genomes
    for enum,file_ in enumerate(sorted(os.listdir(genome_files_dir))):
        # if the user did not supply the name of a reference to use in progressivemauve then make the first file in the directory the reference
        if pmauve_ref == None:
            if enum == 0:
                pmauve_ref = file_
        #/
        pmauve_genome_list.append(genome_files_dir+'/'+file_)
    #/
    # verify that the pmauve reference exists in genomes-directory (only practically effective when the user supplied a name)
    if not os.path.exists(genome_files_dir+'/'+pmauve_ref):
        print('WARNING: could not locate file for reference genome during progressivemauve: '+pmauve_ref)
        print('If you entered this manually, please make sure the file exist in the genome directory: '+str(genome_files_dir))
        print('Terminating!')
        sys.exit()
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
        print('Using '+outgroup_file+' as outgroup',flush=True)
        #/
        # Add outgroup to genome-list
        pmauve_genome_list.append(workdir+'/'+'outgroup'+'/'+outgroup_file)
        #/
    #/
    # Check if number of genomes is below 3, then abort (IQTREE requires at least 3 genomes)
    if len(pmauve_genome_list) < 3:
        print('Unable to generate accurate 3: number of input sequences too low. At least 3 genomes (outgroup included, if applicable) are needed to run IQ-TREE')
        print('Current execution has: '+str(len(pmauve_genome_list))+ ' genomes.')
        sys.exit()
    #/
    # soft-link genomes
    print('Soft-linking genomes...',flush=True)
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
        query_xmfa_out = output_path+'/'+query_name+'.xmfa' #pmaouve output
        query_fa_out = output_path+'/'+query_name+'.fa' #x2fa output
        
        # Run progressiveMauve
        print('\nExecuting progressiveMauve',flush=True)
        pmauve_cmd = ['progressiveMauve','--output',query_xmfa_out,reference_path,query_path]
        subprocess.call(' '.join(map(str,pmauve_cmd)),shell=True)
        #/
        
        # Run x2fa
        print('\nExecuting x2fa.py')
        phylo_file_loc = os.path.abspath(__file__)
        x2fa_file_loc = os.path.dirname(phylo_file_loc) + '/' + 'x2fa.py'
        pmauve_cmd = ['python2',x2fa_file_loc,query_xmfa_out,reference_path,0,query_fa_out]
        subprocess.call(' '.join(map(str,pmauve_cmd)),shell=True)
        #/
        
        # Check if output file was created, else return false
        pmauve_has_output = False
        if os.path.exists(query_fa_out) and os.path.getsize(query_fa_out) > 0:
            pmauve_has_output = True
        #/
        if print_status:            print(print_status + ' finished!',flush=True)
        return pmauve_has_output
    #/
    # setup multiprocessing
    pmauve_outdir = phylo_wd+'/'+'progressivemauve'
    os.makedirs(pmauve_outdir)
    print('Processing genomes (number of parallel processes is '+str(num_threads)+ ' at '+str(1)+' threads each)',flush=True) # progressiveMauve runs single-threaded
    #/
    # start jobs
    if progressivemauve_single:
        print('Argument --progressivemauve_single applied, will not parallelize progressivemauve',flush=True)
        pool = Pool(processes=1)
    else:
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

    print('\tProcessing done!',flush=True)
    
    if not all(list(jobs_status.values())):
        print('Some jobs failed. Terminating.')
        sys.exit()
    #/
    ##/
    ### Concatenate single-alignment-fasta to multi-alignment-fasta and fix fasta-headers
    print('Concatenating alignment-fasta files to MSA-fasta...',flush=True)
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
    
    # Check if ignore IQTREE (--halt_after_msa enabled)
    skip_IQTREE = False
    if halt_after_msa:
        skip_IQTREE = True
    #/
    
    ### IQTREE-related
    if not skip_IQTREE:
        ## Run IQTREE (inside phylo_wd, iqtree produces some file)
        cmd_iqtree = ['iqtree','-T','auto','--threads-max',num_threads,'-s',phylo_wd+'/'+'msa.fasta']
        # check if using outgroup
        if not skip_outgroup:
            cmd_iqtree += ['-o',outgroup_file.replace('.fasta','')]
        #/
        # check if using bootstrap
        if bootstrap_num_iterations > 0:
            cmd_iqtree += ['-B',bootstrap_num_iterations]
        #/
        subprocess.call(' '.join(map(str,cmd_iqtree)),shell=True)
        ##/
        
        ## Save which sample was used as reference in Progressivemauve and copy-out tree-file
        with open(phylo_wd+'/'+'tree_reference.txt','w') as nf:
            nf.write(pmauve_ref+'\n')
        
        with open(phylo_wd+'/'+'msa.fasta.treefile','r') as f:
            with open(output_file_path,'w') as nf:
                for line in f:
                    nf.write(line)
        ##/
    ###/
    ### If user wants to keep tree, copy it to outer
    if keep_msa:
        try:
            shutil.copy2(phylo_wd+'/'+'msa.fasta',msa_output_file_path)
        except Exception as e:
            print('Was unable to copy-out multiple sequence alignment to workdir! (--keep_msa)')
            print(e)
            print('ignoring error and proceeding.')
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
        print('Using '+outgroup_file+' as outgroup',flush=True)
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
