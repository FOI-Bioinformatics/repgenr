#!/usr/bin/env python3

"""
Software name: FlexGenR (Flexible Genome Repository)
Software name: RepGenR (Representative-Genome Repositories)
"""

import sys
import os
import subprocess
import argparse

# Take input arguments
args = sys.argv
if not len(args) >= 2:
    args.append('') # dummy-add something to trigger error-print below
    
wrapper = args[0]
module = args[1]
module_args = args[2:]
if not module in ('metadata','genome','derep','phylo','tree2tax',):
    print('Usage: '+wrapper+' module {metadata,genome,derep,phylo,tree2tax}')
    if module != '':
        print('error: module "'+module.replace('.py','')+'" does not exist!')
    sys.exit()

cmd = [module+'.py'] + module_args

print('[DEBUG] Wrapper CMD:')
print(cmd)
print()

if 0:
    os.system(' '.join(cmd))

if 1:
    subprocess.call(' '.join(cmd),shell=True)

if 0:
    process = subprocess.Popen(' '.join(cmd),shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    while process.poll() is None:
        line_raw = process.stdout.readline() # bufsize?
        line = line_raw.decode()
        line = line.strip('\n')
        if not line: continue
        print(line)
    
"""
process = subprocess.Popen('python download.py -wd output_test2', shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#for line in p.stdout.readlines():
#    line = line.decode()
#    print(line)
while process.poll() is None:
    line_raw = process.stdout.readline() # bufsize?
    line = line_raw.decode()
    line = line.strip('\n')
    if not line: continue

    print(line)
"""



