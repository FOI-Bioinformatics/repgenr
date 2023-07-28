#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import filecmp


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-n','--name',required=False,help='Name to use for packing/unpacking of derep run')
parser.add_argument('--list',action='store_true',help='If specified, will list saved runs')
parser.add_argument('--pack',action='store_true',help='If specified, will save derep files under provided name')
parser.add_argument('--unpack',action='store_true',help='If specified, will load derep files from provided name')
parser.add_argument('--delete',action='store_true',help='If specified, will delete the saved runs')
#/
# parse input
args = parser.parse_args()

workdir = args.workdir

run_name = args.name

list_saved = args.list
execute_pack = args.pack
execute_unpack = args.unpack
delete_saved = args.delete
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory at:')
    print(workdir)
    sys.exit()

actions_given = 0
for i in (list_saved,execute_pack,execute_unpack,delete_saved):
    if i:
        actions_given += 1
if not actions_given == 1:
    print('No action or multiple actions supplied! Please see --help for use instructions and supply one action (e.g., --list or --pack or --unpack)')
    sys.exit()
    
if (execute_pack or execute_unpack or delete_saved) and not run_name:
    print('A run name must be supplied to pack/unpack/delete!')
    sys.exit()
#/
###/

### Define path to stored runs
stored_runs_path = workdir+'/'+'derep_stocker'
###/

### Define expected derep files
expected_derep_files = ['derep_parameters.txt','derep_clustered_genomes.tsv','derep_genomes_status.tsv','genomes_derep_representants']
###/

### Check if user want to list runs
if list_saved:
    # check if folder for runs does not exist
    if not os.path.exists(stored_runs_path):
        print('No runs have been stored!')
        sys.exit()
    #/
    # parse runs stored
    runs_stored = []
    for file_ in os.listdir(stored_runs_path):
        file_path = stored_runs_path+'/'+file_
        if os.path.isdir(file_path):
            runs_stored.append(file_)
    #/
    # print stored runs
    if not runs_stored:
        print('No runs have been stored!')
    else:
        print('Available runs:')
        for i in sorted(runs_stored):
            print(i)
    #/
    # halt execution when "--list" finished
    sys.exit()
    #/
###/

### Define path to input run-name
run_path = stored_runs_path + '/' + run_name
###/

### Check if user want to delete a run
if delete_saved:
    # check if folder for run does not exist
    if not os.path.exists(run_path):
        print('Unable to find specified name:')
        print(run_name)
        sys.exit('Terminated')
    #/
    # delete run
    shutil.rmtree(run_path)
    print('Deletion successful')
    #/
    # halt execution then --delete finished
    sys.exit()
    #/
###/

### Check if user wants to store a run
if execute_pack:
    # Check if a previous save exist with current name, and if so, prompt to delete it
    if os.path.exists(run_path):
        usr_inp = input('WARNING: a previous save exist for this run. Do you want to overwrite it? (y/n) ')
        if usr_inp.lower() in ('y','yes',):
            shutil.rmtree(run_path)
        else:
            sys.exit('\nTerminated')
    #/
    # Check if there are drep files in --workdir
    expected_files_found = []
    for expected_file in expected_derep_files:
        if os.path.exists(workdir+'/'+expected_file):
            expected_files_found.append(True)
        else:
            expected_files_found.append(False)
    if not all(expected_files_found):
        print('FATAL: did not find the expected files from derep. Please make sure you have executed the derep module.')
        print('Expected files:')
        for i in expected_derep_files:
            print(i)
        sys.exit('\nTerminated')
    #/
    # Init dirs
    if not os.path.exists(stored_runs_path):            os.makedirs(stored_runs_path)
    if not os.path.exists(run_path):                    os.makedirs(run_path)
    #/
    ## Save derep files (excluding dereplicated genomes; will soft-link in those)
    # save files
    for expected_file in expected_derep_files:
        if os.path.isfile(workdir+'/'+expected_file):
            shutil.copy2(workdir+'/'+expected_file,run_path+'/'+expected_file)
    #/
    # link-in genomes listed in {--workdir}/genomes_derep_representants from {--workdir}/genomes into {--workdir}/derep_stocker/{run_name}/genomes_derep_representants
    os.makedirs(run_path+'/'+'genomes_derep_representants')
    for file_ in os.listdir(workdir+'/'+'genomes_derep_representants'):
        os.symlink(os.getcwd()+'/'+workdir+'/'+'genomes'+'/'+file_,run_path+'/'+'genomes_derep_representants'+'/'+file_)
    #/
    ##/
    # halt execution when --pack finished
    print('Packing of run completed!')
    sys.exit()
    #/
###/

### Check if user wants to load a run
if execute_unpack:
    # check if input run exists and have the expected files
    if not os.path.exists(run_path):
        print('Unable to find specified name:')
        print(run_name)
        sys.exit('Terminated')
    
    expected_files_found = []
    for expected_file in expected_derep_files:
        if os.path.exists(run_path+'/'+expected_file):
            expected_files_found.append(True)
        else:
            expected_files_found.append(False)
    if not all(expected_files_found):
        print('FATAL: did not find the expected files from derep to load. Please make sure you have provided the correct name of a previously saved derep run')
        print('Expected files:')
        for i in expected_derep_files:
            print(i)
        sys.exit('\nTerminated')
    #/
    ## Check if the current run under {--workdir} is saved (we will load a saved run to replace the current one)
    # Check if expected files exist in --workdir
    expected_files_found = []
    for expected_file in expected_derep_files:
        if os.path.exists(workdir+'/'+expected_file):
            expected_files_found.append(True)
        else:
            expected_files_found.append(False)
    #/
    # If they exist, then check if they have been saved and prompt to halt execution
    if any(expected_files_found) and not all(expected_files_found):
        print('WARNING: derep files under --workdir appears to be corrupted. Can not check if they were saved previously. Please retry after manually removing these files:')
        for i in sorted(expected_derep_files):
            print(i)
        sys.exit('\nTerminated')
    
    if all(expected_files_found):
        existing_saves_file_match = []
        for existing_save in os.listdir(stored_runs_path):
            expected_files_are_identical = []
            existing_save_path = stored_runs_path + '/' + existing_save
            for expected_file in expected_derep_files:
                # Skip directory (comparing folder of files [in --workdir] to folder of symlinks [run_path]. will implement a check later if needed)
                if os.path.isdir(workdir+'/'+expected_file): continue
                #/
                # Check if files are identical
                files_are_identical = filecmp.cmp(workdir+'/'+expected_file, existing_save_path+'/'+expected_file,shallow=True)
                expected_files_are_identical.append(files_are_identical)
                #/
            
            # Check if all files were identical
            if all(expected_files_are_identical):
                existing_saves_file_match.append(existing_save)
            #/
        
        # check if we did not have a save of current run
        if not existing_saves_file_match:
            usr_inp = input('WARNING: did not have a save of current files in --workdir folder. Do you want to proceed loading without saving? This will overwite the current derep-files in --workdir. (y/n) ')
            if not usr_inp.lower() in ('y','yes',):
                sys.exit('\nTerminated')
        #/
    #/
    ## Load run
    # copy-out files
    for expected_file in expected_derep_files:
        if os.path.isfile(run_path+'/'+expected_file):
            shutil.copy2(run_path+'/'+expected_file,workdir+'/'+expected_file)
    #/
    # copy representative genomes in {run_path}/genomes_derep_representants to {--workdir}/genomes_derep_representants from {--workdir}/genomes (they are symlinks in {run_path}/genomes_derep_representants)
    #@ clean previous
    if os.path.exists(workdir+'/'+'genomes_derep_representants'):       shutil.rmtree(workdir+'/'+'genomes_derep_representants')
    #@/
    #@ init dir
    os.makedirs(workdir+'/'+'genomes_derep_representants')
    #@/
    #@ populate dir
    for file_ in os.listdir(run_path+'/'+'genomes_derep_representants'):
        shutil.copy2(workdir+'/'+'genomes'+'/'+file_,workdir+'/'+'genomes_derep_representants'+'/'+file_)
    #@/
    #/
    # halt execution when --unpack finished
    print('Unpacking of run completed!')
    sys.exit()
    #/
###/
