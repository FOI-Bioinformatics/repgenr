#!/usr/bin/env python3

import sys
import os
import subprocess
import argparse
import time


submodules_available = ('metadata','vmetadata','genome','vgenome','glance','derep','derep_unpack','derep_stocker','phylo','tree2tax',)
software_description = '''
                        RepGenR:
                            
                        '''

## Define argparse
argparser = argparse.ArgumentParser(description=software_description)
argparser.add_argument('submodule', choices=submodules_available, help='Submodule to run')
argparser.add_argument('args', nargs=argparse.REMAINDER, help='Arguments to pass to the submodule')
##/
## Check if had any user input
if len(sys.argv) == 1:
    argparser.print_help(sys.stderr)
    sys.exit()
##/
## Parse user input
args = argparser.parse_args()
submodule = args.submodule
submodule_args = args.args
##/

## Validate input
# Check if arguments are input multiple times
args_count = {}
for arg in submodule_args:
    if not arg in args_count:       args_count[arg] = 0
    args_count[arg] += 1

args_count_multiple = []
for arg,count in args_count.items():
    if arg[0] == '-' or arg[:2] == '--':
        if count >= 2:
            args_count_multiple.append([arg,count])

if args_count_multiple:
    print('[WARNING] Found multiple inputs of arguments:')
    for arg,count in args_count_multiple:
        print(arg+'\t'+str(count))
    print(' '.join(args))
    print('Please specify arguments once. Terminating!')
    sys.exit()
#/
##/
### Execute submodule
## Construct the command to execute the submodule
# init cmd with submodule call
cmd = [submodule + '.py']
#/
# append arguments to the submodule. Check if the arguments have spaces within then, then add quotation
submodule_args2 = []
for entry in submodule_args:
    # add qutoes around inputs with spaces
    if entry.find(' ') != -1:
        entry = '"'+entry+'"'
    #/
    # save
    submodule_args2.append(entry)
    #/
cmd += submodule_args2
#/
# Construct logging module: merge stdout and stderr (captures output of all subprocesses and subcommands) and save to log through "tee". Try to get the workdir passed to subcommand
workdir = ''
for arg_enum,_ in enumerate(submodule_args):
    if arg_enum == 0: continue # backawards checking
    argument = submodule_args[arg_enum-1]
    value = submodule_args[arg_enum]
    
    if argument in ('--workdir','-wd',):
        if submodule == 'metadata' or os.path.exists(value): # dont mess upp logging if a non-existant directory is passed. If the subcommand is "metadata", first subcommand to call, there will be no directory for -wd.
            workdir = value
            break # stop on first

log_cmd_extension = []
if workdir:
    if submodule == 'metadata' and not os.path.exists(workdir):     os.makedirs(workdir) # if subcommand is metadata, then make the directory for log
    log_cmd_extension = ['2>&1 | tee -a '+workdir+'/'+'repgenr.log']

cmd += log_cmd_extension
#/
##/
# Write command passed to log
if log_cmd_extension:
    with open(workdir+'/'+'repgenr.log','a') as fo:
        fo.write('[CMD] '+' '.join(cmd)+'\n')
        fo.write('[CWD] '+os.getcwd()+'\n')
#/
# Run the subprocess with the constructed command
cmd_timestamp_start = ['printf "[INFO][TIMESTAMP][START] %s\n" "$(date)"']
cmd_timestamp_start += log_cmd_extension
subprocess.call(' '.join(cmd_timestamp_start),shell=True)

time_start = time.time()
subprocess.call(' '.join(cmd),shell=True)
time_end = time.time()

cmd_timestamp_end = ['printf "[INFO][TIMESTAMP][END] %s\n" "$(date)"']
cmd_timestamp_end += log_cmd_extension
subprocess.call(' '.join(cmd_timestamp_end),shell=True)

time_delta = time_end - time_start
mins,secs = divmod(time_delta,60)
hours,mins = divmod(mins,60)
cmd_runtime = ['printf "[INFO][RUNTIME] '+f'Elapsed time: {int(hours):02d} hours and {int(mins):02d} minutes\n"']
cmd_runtime += log_cmd_extension
subprocess.call(' '.join(cmd_runtime),shell=True)
#/
###/
# add dummy function for setup.py entry-point
def main():
    pass
#/
