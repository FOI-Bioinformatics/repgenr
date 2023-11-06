#!/usr/bin/env python3

import os
import sys
import argparse
import zipfile
import gzip
import ast
import shutil
import subprocess
import time


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata-command')
parser.add_argument('--accession_list_only',action='store_true',help='If specified, will output NCBI download accession list and then terminate')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
#/
# parse input
args = parser.parse_args()

workdir = args.workdir
halt_after_accession_list = args.accession_list_only
keep_files = args.keep_files
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (created by metadata-command). Please check the input:')
    print(workdir)
    sys.exit()
if not os.path.exists(workdir+'/'+'outgroup_accession.txt'):
    print('Could not locate file outgroup_accession.txt to use for tree rooting (created by metadata-command)')
    sys.exit()
#/
###/

### Parse metadata
## selected accessions metadata
try:
    print('Parsing metadata...')
    with open(workdir+'/'+'metadata_selected.tsv','r') as f:
        line = f.readline()
        accessions = ast.literal_eval(line)
except:
    print('Could not read metadata at location:')
    print(workdir+'/'+'metadata_selected.tsv')
    print('Have you run the metadata-command?')
    sys.exit()
##/
## Outgroup metadata
try:
    with open(workdir+'/'+'metadata_outgroup.tsv','r') as f:
        line = f.readline()
        outgroup_metadata = ast.literal_eval(line)
except:
    print('Could not read metadata at location:')
    print(workdir+'/'+'metadata_outgroup.tsv')
    print('Have you run the metadata-command?')
    sys.exit()
##/
###/



### Download accessions
def format_acc_output_file(acc_data_entry):
    """
    Input: a metadata entry, will parse the 'tax_gtdb' field and the 'accession' field
    Output: returns name on format {family}_{genus}_{species}_{GCF/GCA}.fasta
    """
    tax = acc_data_entry['tax_gtdb']
    return tax['family']+'_'+tax['genus']+'_'+tax['species']+'_'+acc_data_entry['accession']+'.fasta'

## Remove accessions that do not exist in the input
if os.path.exists(workdir+'/'+'genomes'):
    for file_ in os.listdir(workdir+'/'+'genomes'):
        fam,gen,spec,acn_type,acn_plusExtension = file_.split('_')
        acn = acn_plusExtension.replace('.fasta.gz','').replace('.fasta','')
        
        accession = acn_type+'_'+acn
        
        if not accession in accessions:
            os.remove(workdir+'/'+'genomes'+'/'+file_)
            print('Removed '+ file_ +' (no longer exists in input)')
##/

## Get accessions that have not been downloaded yet
accessions_to_download = []
for acc,data in accessions.items():
    # check if we already downloaded this genome accession (and that file is not empty)
    acc_genome_file = format_acc_output_file(data)
    if not os.path.exists(workdir+'/'+'genomes'+'/'+acc_genome_file) or os.path.getsize(workdir+'/'+'genomes'+'/'+acc_genome_file) == 0:
        accessions_to_download.append(acc)
    #/
print('Will try to download '+str(len(accessions_to_download))+' genomes...')
##/

## Dump list of accessions for "ncbi datasets" software
with open(workdir+'/'+'ncbi_acc_download_list.txt','w') as nf:
    nf.write('\n'.join(accessions_to_download))

if halt_after_accession_list:
    print('Accession list printed. Terminating')
    sys.exit()
##/

## Download genome fastas
def ncbi_downloader(ncbi_datasets_cmd):
    #
    # This function passes an input command to NCBI datasets through Popen and parses its outout in a "clean" way
    #
    
    prev_write_timestamp = None
    rows_100_perc_written = set()
    try:
        # print execution start
        print('[NCBI-DATASETS] started',flush=True)
        #/
        ncbi_datasets_process = subprocess.Popen(ncbi_datasets_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        while ncbi_datasets_process.poll() is None:
            line_raw = ncbi_datasets_process.stdout.readline() # bufsize?
            line = line_raw.decode()
            line = line.strip('\n')
            if not line: continue
            
            # write status update lines
            if (prev_write_timestamp == None) or (time.time() - prev_write_timestamp >= 5) or (line.find('100%') != -1): #only print once every X:th second (log gets clogged otherwise) but allow "finish/100%"-line
                # check if "finish/100%"-line and if it was already printed
                if line.find('100%') != -1:
                    # try to remove any ANSI-code that tells the cursor to move up "[2K"
                    try:
                        if line.find('[2K') != -1:
                            line=line.split('[2K')[1]
                    except:
                        pass
                    #/
                    if line in rows_100_perc_written: continue # skip if already printed
                    rows_100_perc_written.add(line) # add to memory so it is not printed again
                #/
                prev_write_timestamp = time.time()
                if line.startswith('Downloading') or (line.find('Collecting') != -1 and line.find('records') != -1):
                    print('[NCBI-DATASETS] ' + line,flush=True)
                    pass # I was going to do something else here. For now, all lines are printed
                else:
                    print('[NCBI-DATASETS] ' + line,flush=True)
            #/
        # print final/flush
        print('[NCBI-DATASETS] finished',flush=True)
        #/
    except:
        print('Failed to download from NCBI using "datasets"-software. You can try to run the command manually:')
        print(' '.join(ncbi_datasets_cmd))
        sys.exit()

print('Running NCBI-DATASETS software (download genomes)')
if accessions_to_download:
    ncbi_datasets_cmd = ['datasets','download','genome','accession',
                         '--inputfile',workdir+'/'+'ncbi_acc_download_list.txt',
                         '--filename',workdir+'/'+'ncbi_download.zip']
    
    ncbi_downloader(ncbi_datasets_cmd)
    
    ### Unpack genomes from downloaded zip (into formatted gzipped files)
    print('Unpacking genomes from NCBI download into "genomes"-folder...')
    try:
        with zipfile.ZipFile(workdir+'/'+'ncbi_download.zip','r') as zip_fo:
            if not os.path.exists(workdir+'/'+'genomes'):        os.makedirs(workdir+'/'+'genomes')
            try:
                zip_fo.extractall(workdir+'/'+'ncbi_extract')
            except:
                print('ZIP extractall failed. Possibly, it is due to a corrupted download from NCBI. Try re-running or else raise an issue. Terminating!')
                sys.exit()
            
            # subdir will be at ../ncbi_extract/ncbi_dataset/data/GCx_XXXXXX in extracted folder
            for path,dirnames,filenames in os.walk(workdir+'/'+'ncbi_extract'):
                for dirname in dirnames:
                    if dirname in accessions:
                        # Get file from directory
                        genome_fa_path = None
                        for file_ in os.listdir(path+'/'+dirname):
                            if file_.endswith('.fna'):
                                genome_fa_path = path+'/'+dirname+'/'+file_
                                break
                        
                        if not genome_fa_path:
                            print('Could not locate genome fasta for downloaded accession: '+dirname)
                        
                        # copy genome file to genomes folder
                        acc_genome_file = format_acc_output_file(accessions[dirname])
                        genome_fa_target_path = workdir+'/'+'genomes'+'/'+acc_genome_file
                        shutil.copy(genome_fa_path,genome_fa_target_path)
                        #/
    except zipfile.BadZipFile:
        print('ZIP-file from NCBI download appears to be corrupted - was the connection interrupted? Try re-running this module to see if network issues persist.')
        sys.exit()
    except FileNotFoundError:
        print('ZIP-file from NCBI download was not found - was the connection interrupted? Try re-running this module to see if network issues persist.')
        sys.exit()
            
    # remove ZIP file
    if not keep_files:
        print('Cleaning workspace...')
        os.remove(workdir+'/'+'ncbi_download.zip')
        shutil.rmtree(workdir+'/'+'ncbi_extract')
        #os.remove(workdir+'/'+'ncbi_acc_download_list.txt')
    #/
else:
    print('All genomes were already downloaded!')
###/

### Download outgroup sample
print('Downloading and extracting outgroup sample...')
outgroup_accession = None
with open(workdir+'/'+'outgroup_accession.txt','r') as f:
    outgroup_accession = f.readline().strip('\n')

ncbi_datasets_cmd_outgroup = ['datasets','download','genome','accession',
                              outgroup_accession,
                              '--filename',workdir+'/'+'ncbi_download_outgroup.zip']

ncbi_downloader(ncbi_datasets_cmd_outgroup)

zip_fo = zipfile.ZipFile(workdir+'/'+'ncbi_download_outgroup.zip','r')
if not os.path.exists(workdir+'/'+'outgroup'):        os.makedirs(workdir+'/'+'outgroup')
for item in zip_fo.namelist():
    if item.endswith('.fna'): # in NCBI zip-file, the genomes ends with this
        dirname = item.split('/')[-2]
        basename = item.split('/')[-1]
        
        if dirname == outgroup_accession:
            acc_genome_file = format_acc_output_file(outgroup_metadata[dirname])
            genome_fa_target_path = workdir+'/'+'outgroup'+'/'+acc_genome_file
            with open(genome_fa_target_path,'w') as nf:
                nf.write(zip_fo.read(item).decode('utf-8'))
            #/
zip_fo.close()

os.remove(workdir+'/'+'ncbi_download_outgroup.zip')
print('Cleaning workspace...')
###/