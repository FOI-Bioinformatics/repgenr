#!python

import sys
import os
import subprocess
import argparse
import time


submodules_available = ('metadata','genome','glance','derep','derep_stocker','phylo','tree2tax',)
software_description = 'RepGenR:'

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
## Execute submodule
# Construct the command to execute the submodule
cmd = [submodule + '.py'] + submodule_args
#/
# Construct logging module: merge stdout and stderr and save to log through "tee". Try to get the workdir passed to subcommand
workdir = ''
for arg_enum,_ in enumerate(submodule_args):
    if arg_enum == 0: continue # backawards checking
    argument = submodule_args[arg_enum-1]
    value = submodule_args[arg_enum]
    
    if argument in ('--workdir','-wd',):
        if os.path.exists(value):
            workdir = value
            break # stop on first

if workdir:
    log_cmd_extension = ['2>&1 | tee -a '+workdir+'/'+'repgenr.log']
else:
    log_cmd_extension = []
    
cmd += log_cmd_extension
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
##/
# add dummy function for setup.py entry-point
def main():
    pass
#/
