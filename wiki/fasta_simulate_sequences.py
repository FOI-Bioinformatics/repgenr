#!/data/users/jaclew/software/anaconda3/bin/python

import argparse
import gzip
import sys
from random import randint

###### This script will import a reference-fasta file and then output randomized strings from that reference. #######

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='Produces random strings (fasta) from an input sequence (fasta) to stdout\nSee --help for options')
parser.add_argument('input',help='Path to working directory, created by metadata-command')
parser.add_argument('-l','--out_length',type=int,default=1000,help='Sets the attempted length of output sequences (default: 1000)')
parser.add_argument('-n','--num_seqs',type=int,default=1000,help='Sets the number of output sequences (default: 1000)')
parser.add_argument('-b','--basename',help='Sets the basename for read names (e.g. --basename read outputs: read1, ..., readN)')
#/
# parse input
args = parser.parse_args()

input_fasta = args.input
out_length = args.out_length
num_seqs = args.num_seqs
read_basename = args.basename
#/
###/

## Parse sequences
file_is_gzipped = False
if input_fasta.endswith('.gz'):
    f = gzip.open(input_fasta,'rb')
    file_is_gzipped = True
else:
    f = open(input_fasta,'r')
    
#import sequences
rn = f.readline()
if file_is_gzipped: rn = rn.decode() #handle if zipped

if rn[0] == '>':
    fastx = 'fasta'
elif rn[0] == '@':
    fastx = 'fastq'
else:
    sys.exit('this is not a fasta? [first entry does not have > at line start]')
rn = rn.strip('\n')[1:]

ref_seqs = {}
seq = ''
for line in f:
    if file_is_gzipped: line = line.decode() #handle if zipped
    
    line = line.strip('\n')
    
    #keep storing seq to current RN
    if not line[0] == '>':
        seq += line
    #if we find new rn, store old and write it down along with folded seq
    else:
        line = line[1:]
        rn_old = rn
        rn = line
        
        #handle sequence
        ref_seqs[rn_old] = seq
        seq = '' #reset
        

#/import  sequences
#import last
if seq:
    rn_old = rn
    ref_seqs[rn_old] = seq
    seq = '' #reset
#/import last
##/

## Print random seq
refs = list(ref_seqs.keys())
output_seqs = 0
while output_seqs < num_seqs:
    
    rand_ref = refs[randint(0,len(refs)-1)]
    ref_seq = ref_seqs[rand_ref]
    
    rand_pos = randint(0,len(ref_seq))
    rand_seq = ref_seq[rand_pos:rand_pos+out_length]
    
    if not rand_seq: continue # skip if no seq was parsed (should not happen)
    
    out_header = rand_ref.split()[0]
    if read_basename:       out_header = read_basename
    out_header = out_header+'_'+str(output_seqs)
    
    print('>'+out_header+'\n'+rand_seq)
    output_seqs += 1
##/