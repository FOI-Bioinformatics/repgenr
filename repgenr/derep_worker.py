#!/usr/bin/env python3

### IDE: this script used to be the "derep.py". The new "derep.py" is a wrapper around this script that enable chunking of large datasets.

import argparse
import os
import sys
import subprocess
import shutil


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-wd','--workdir',required=True,help='Path to working directory, created by metadata and genome commands')
parser.add_argument('-sani','--secondary_ani',default=0.99,help='Average nucleotide identity (ANI) threshold for clustering in sensitive step of dRep (secondary clustering)')
parser.add_argument('-pani','--primary_ani',default=0.90,help='Average nucleotide identity (ANI) threshold for clustering in rough step of dRep (primary clustering)')
parser.add_argument('--S_algorithm',default='fastANI',help='Algorithm to use in dRep fastANI. Possible values as of dRep v3.4.1: fastANI,gANI,goANI,ANIn,ANImf (default: fastANI)')
parser.add_argument('-t','--threads',type=int,default=24,help='Number of threads to use')
parser.add_argument('-l','--length',type=int,default=0,help='Length threshold for sequences')
parser.add_argument('-nc','--cov_thresh',type=float,default=None,help='See dRep documentation')
parser.add_argument('-N50W','--N50_weight',type=float,default=None,help='See dRep documentation')
parser.add_argument('-sizeW','--size_weight',type=float,default=None,help='See dRep documentation')
parser.add_argument('--ignoreGenomeQuality',action='store_true',default=None,help='See dRep documentation')
parser.add_argument('--clusterAlg',type=str,default=None,help='See dRep documentation')
#/
# parse input
args = parser.parse_args()

workdir = args.workdir

secondary_ani = args.secondary_ani
primary_ani = args.primary_ani
S_algorithm = args.S_algorithm
num_threads = args.threads

cov_thresh = args.cov_thresh
N50_weight = args.N50_weight
size_weight = args.size_weight
ignoreGenomeQuality = args.ignoreGenomeQuality
clusterAlg = args.clusterAlg
seq_len_threshold = args.length
#/
# validate input
if not os.path.exists(workdir):
    print('Could not locate working directory (make sure you have run metadata and genome commands):')
    print(workdir)
    sys.exit()
#/
###/

### Must decompress all gzipped files (for prodigal, that is run within DREP)
print('Decompressing genomes (required for prodigal, that is run within DREP)...',flush=True)
# Check if there are files to decompress..
run_decompress = False
for file_ in os.listdir(workdir+'/genomes'):
    if file_.endswith('.gz'):
        run_decompress = True
        break
#/
# Run decompress
if run_decompress:
    cmd_decompress = ['gzip','-d',workdir+'/genomes/*.gz']
    subprocess.call(' '.join(map(str,cmd_decompress)),shell=True)
#/
###/

### Run dRep (sa = ANI threshold to form secondary clusters [ALIGNMENT]. pa = ANI thresholds to form primary clusters [MASH])
print('Running genome dereplication using dRep...',flush=True)
cmd_drep = ['dRep','dereplicate',workdir+'/drep_workdir','-g',workdir+'/genomes/*.fasta','--processors',num_threads,'-sa',secondary_ani,'-pa',primary_ani,'--S_algorithm',S_algorithm,'--length',seq_len_threshold]

# check if append virus-related arguments
if cov_thresh != None:             cmd_drep += ['--cov_thresh',cov_thresh]
if N50_weight != None:             cmd_drep += ['--N50_weight',N50_weight]
if size_weight != None:            cmd_drep += ['--size_weight',size_weight]
if ignoreGenomeQuality != None:    cmd_drep += ['--ignoreGenomeQuality']
if clusterAlg != None:             cmd_drep += ['--clusterAlg',clusterAlg]
#/

cmd_drep.append('-d') #debug output
subprocess.call(' '.join(map(str,cmd_drep)),shell=True)

if not os.path.exists(workdir+'/drep_workdir'):
    print('Could not locate dRep working directory! Please confirm that dRep is installed and can execute properly.')
    sys.exit()
###/

### Copy repesentative genomes
print('Copying representing genomes from dereplication...',flush=True)
if not os.path.exists(workdir+'/genomes_derep_representants'):       os.makedirs(workdir+'/genomes_derep_representants')
cmd_cp = ['cp','-r',workdir+'/drep_workdir/dereplicated_genomes/*',workdir+'/genomes_derep_representants']
subprocess.call(' '.join(map(str,cmd_cp)),shell=True)
###/

### Get info about representative sequences and contained sequences and failed genomes
genomes_status = {} # genome -> "status in [representative,contained,failed]"

## Parse genome representatives
genome_representatives = {} # representative_file_name -> cluster (cluster is parsed below)
for file_ in os.listdir(workdir+'/'+'drep_workdir/dereplicated_genomes'):
    genome_representatives[file_] = None
##/
## Parse cluster-contained genomes
clusters_genomes = {}
with open(workdir+'/'+'drep_workdir/data_tables/Cdb.csv','r') as f:
    for enum,line in enumerate(f):
        if enum == 0: continue #skip header
        line = line.strip('\n')
        line = line.split(',')
        genome,clust = line[0],line[1]
        
        if not clust in clusters_genomes:        clusters_genomes[clust] = set()
        clusters_genomes[clust].add(genome)
        
        # check if current line is the representative
        if genome in genome_representatives:
            genome_representatives[genome] = clust
            genomes_status[genome] = 'representative'
        else:
            genomes_status[genome] = 'contained'
        #/
##/
## Check if we run with virus-settings, then modify genomeInformation-file.
## for viral run, columns ARE REPORTED DIFFERENTLY as of using --ignoreGenomeQuality
## It reports columns "genome_path,_length,_N50,genome,centrality"
# Check if there are 5 columns only in file
genomeInformation_is_viral = False
with open(workdir+'/'+'drep_workdir/data_tables/genomeInformation.csv','r') as f:
    for enum,line in enumerate(f):
        line = line.strip('\n')
        line = line.split(',')
        
        if len(line) != 7:
            genomeInformation_is_viral = True
#/
# Copy-away genomeInformation file and write a new one on same format as "default" dRep
if genomeInformation_is_viral:
    print('INFO: genomeInformation identified as 5-column format. Will make backup and re-create it as the expected 7-column file',flush=True)
    shutil.move(workdir+'/'+'drep_workdir/data_tables/genomeInformation.csv',workdir+'/'+'drep_workdir/data_tables/genomeInformation_orig.csv')
    
    with open(workdir+'/'+'drep_workdir/data_tables/genomeInformation_orig.csv','r') as f:
        with open(workdir+'/'+'drep_workdir/data_tables/genomeInformation.csv','w') as nf:
            for enum,line in enumerate(f):
                # check if header. Then write the 7-column header instead
                if enum == 0:
                    writeArr = ['genome','completeness','contamination','strain_heterogeneity','length','N50','centrality']
                    nf.write(','.join(writeArr)+'\n')
                    continue
                #/
                # parse original line and information
                line = line.strip('\n')
                line = line.split(',')
                genome_path,length,N50,genome,centrality = line
                #/
                # write modofied line with 7-column format structured by the expected columns
                completeness = 'NA'
                contamination = 'NA'
                strainHeterogeneity = 'NA'
                writeArr = [genome,completeness,contamination,strainHeterogeneity,length,N50,centrality]
                nf.write(','.join(map(str,writeArr))+'\n')
                #/
##/
## Get genomes that was not involved in dereplication (i.e., failed checkM or similar) NOTE: these genomes can still be representatives if no better genome existed
with open(workdir+'/'+'drep_workdir/data_tables/genomeInformation.csv','r') as f:
    for enum,line in enumerate(f):
        # Expected columns:
        # genome	completeness	contamination	strain_heterogeneity	length	N50	centrality
        # centrality appears so be set to 0 for genomes that were not involved in formation of derep clusters
        if enum == 0: continue # skip header
        line = line.strip('\n')
        line = line.split(',')
        genome,_completeness,_contamination,_strainHeterogeneity,_length,_N50,centrality = line
        try:
            if float(centrality) == 0:
                genomes_status[genome] = 'fail_qc'
        except:
            print('Unexpected data in dRep file "genomeInformation". Please open issue.')
            sys.exit()
##/
## Write cluster-contained genomes
with open(workdir+'/'+'derep_clustered_genomes.tsv','w') as nf:
    for genome_rep,clust in genome_representatives.items():
        for genome_contained in clusters_genomes[clust]:
            if genome_contained == genome_rep: continue # skip self
            writeArr = [genome_rep,genome_contained]
            nf.write('\t'.join(writeArr)+'\n')
##/
## Write genome status
with open(workdir+'/'+'derep_genome_status.tsv','w') as nf:
    for genome,status in genomes_status.items():
        writeArr = [genome,status]
        nf.write('\t'.join(writeArr)+'\n')
##/
###/
### Write ANI-thresh used
with open(workdir+'/'+'derep_ANI_threshold.txt','w') as nf:
    nf.write(str(secondary_ani)+'\n')
###/