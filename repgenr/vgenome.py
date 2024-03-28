#!/usr/bin/env python3

import os
import sys
import argparse
import pickle
from Bio import SeqIO
from statistics import mean,median,stdev
import shutil
import subprocess


### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = 
'''
Organize genomes and select outgroup from "vmetadata" module. Virus-equivalent for "genome" module
''')

parser.add_argument('-wd','--workdir',required=True,help='Path to working directory. This folder will be created if not already present')
#parser.add_argument('-t','--taxids',required=False,default=None,help='[NOT EFFECTIVE ATM: select using "Taxnonomy selection" with NCBI] If supplied, will parse sequences by taxonomic identifiers. Multiple values are separated by comma')
#parser.add_argument('-n','--names',required=False,default=None,help='[NOT EFFECTIVE ATM: select using "Taxnonomy selection" with NCBI] If supplied, will parse sequences by taxonomic names (case-independent exact match). Multiple values are separated by comma')
parser.add_argument('--no_outgroup',required=False,action='store_true',default=None,help='If specified, will not determine an outgroup (default: determine outgroup)')
parser.add_argument('--outgroup_candidates_taxid_min_genomes',required=False,type=int,default=5,help='When selecting outgroup candidates, only include taxids with at least this number of genomes (default:5)')


grp_genome_selection = parser.add_argument_group('Genome selection')
grp_genome_selection.add_argument('-f','--filter',required=False,default='complete genome',help='Filter-tag for fasta headers to use (default: "complete genome")')
grp_genome_selection.add_argument('--use_filter',default=None,action='store_true',required=False,help='If specified, will only parse sequences with tag in header (default: do not require tag)')
grp_genome_selection.add_argument('--length_all',default=None,action='store_true',required=False,help='If specified, will select all sequences within the min-max range of genome sequences')
grp_genome_selection.add_argument('--length_deviation',type=int,default=10,required=False,help='Set the percent deviation from median length allowed to keep a sequence (default:10)')
grp_genome_selection.add_argument('--length_method',choices=('median_of_medians','mean_of_means'),default='median_of_medians',required=False,help='Select statistic to use when determining length selection, i.e. "median_of_medians" uses the median sequence length reported per each taxon and takes the median of all taxons (default:"median_of_medians")')
grp_genome_selection.add_argument('--length_range',default=None,required=False,help='Specify the range to select genomes within, e.g. --length 25000-35000 to select genomes that are between 25-35 kbp')
grp_genome_selection.add_argument('--discard',required=False,default=None,help='Discard sequences that match input in the fasta-header. Case-sensitive. Example: "UNVERIFIED,PREDICTED,133337.123" to discard unverified or predicted sequences and sequence with BVBRC identifier "133337.123". Multiple values are separated by comma')

grp_tax_selection = parser.add_argument_group('Taxonomy selection (requires NCBI metadata)')
grp_tax_selection.add_argument('-tg','--target_genus',required=False,default=None,help='Target genus. Strings are interpreted as taxonomy names and integers are interpreted as taxonomy identifiers. Multiple values are separated by comma')
grp_tax_selection.add_argument('-ts','--target_species',required=False,default=None,help='Target species. Strings are interpreted as taxonomy names and integers are interpreted as taxonomy identifiers. Multiple values are separated by comma')
grp_tax_selection.add_argument('-tse','--target_serotype',required=False,default=None,help='Target serotype. Strings are interpreted as taxonomy names and integers are interpreted as taxonomy identifiers. Multiple values are separated by comma')
grp_tax_selection.add_argument('-tc','--target_custom',required=False,default=None,help='Target "custom". Supply a custom taxonomic level and value to match. Example: --target_custom undefined_strain:FX900 to get datasets with FX900 in column "undefined_strain". Multiple values are separated by comma')


#parser.add_argument('--no_ncbi',action='store_true',help='[NOT SUPPORTED ATM: NCBI REQUIRED] If specified, will not use NCBI metadata from previous step (default: use NCBI metadata)')
#parser.add_argument('--limit',type=int,default=None,help='[NOT EFFECTIVE] Limits the analysis to the specified number of genomes (e.g. for test/debug purposes)')
parser.add_argument('--glance',action='store_true',help='If specified, will print dataset selection and terminate')
parser.add_argument('--print_fasta_headers',action='store_true',help='If specified, will print fasta headers of selected sequences. Couple with --glance to terminate before writing genome sequences as output')
parser.add_argument('--keep_files',action='store_true',help='If specified, will save intermediary files')
parser.add_argument('--ignore_duplicates',action='store_true',help='If specified, will prevent termination of the software when faced with multiple identical fasta headers (default: terminate)')
#/
# parse input
if 0 and 'run':
    args = parser.parse_args()
else:
    print('IDE MODE')
    args = parser.parse_args(['--workdir','../../2024/tmp_carro_repgenr/Asfarviridae',
                              '--outgroup_candidates_taxid_min_genomes','0',
                              '--ignore_duplicates',
                              '--target_species','African swine fever virus',
                              '--keep_files'])


workdir = args.workdir

input_taxids_raw = None #input_taxids_raw = args.taxids
input_taxnames_raw = None #input_taxnames_raw = args.names
do_not_generate_outgroup = args.no_outgroup
outgroup_candidates_taxid_min_genomes = args.outgroup_candidates_taxid_min_genomes

# genome selection
tag_required_str = args.filter
tag_required_bool = args.use_filter
seq_len_select_all = args.length_all
seq_len_deviation = args.length_deviation
seq_len_method = args.length_method
seq_len_range = args.length_range
seq_discard_tags = args.discard
#/
# taxonomy selection (requires NCBI metadata)
tax_target_genus = args.target_genus
tax_target_species = args.target_species
tax_target_serotype = args.target_serotype
tax_target_custom = args.target_custom
#/


limit_samples = None # limit_samples = args.limit
no_ncbi = False # no_ncbi = args.no_ncbi
perform_glance = args.glance
perform_fasta_header_print = args.print_fasta_headers
keep_files = args.keep_files
fasta_duplicates_do_not_terminate = args.ignore_duplicates # lets the user disable termination when multiple sequences with the same fasta ID are found (example: >MZ566623)
#/
## Construct input
# "base" input taxids/taxnames
input_taxids = None
input_taxnames = None
if input_taxids_raw:
    input_taxids = [x.strip(' ') for x in input_taxids_raw.split(',')]
if input_taxnames:
    input_taxnames = []
    for entry in input_taxnames.split(','):
        # trim spaces before/after text-content (i.e. if user enters "my species A, my species B" instead of "my species A,my species B")
        input_taxnames.append(entry.rstrip(' ').lstrip(' '))
        #/
#/
# "NCBI metadata" inputs
if tax_target_genus:    tax_target_genus = [x.rstrip(' ').strip(' ').lower() for x in tax_target_genus.split(',')]
if tax_target_species:  tax_target_species = [x.rstrip(' ').strip(' ').lower() for x in tax_target_species.split(',')]
if tax_target_serotype: tax_target_serotype = [x.rstrip(' ').strip(' ').lower() for x in tax_target_serotype.split(',')]
if tax_target_custom:   tax_target_custom = [x.rstrip(' ').strip(' ').lower() for x in tax_target_custom.split(',')]
#/
# Output discard tags (case sensitive)
if seq_discard_tags:    seq_discard_tags = [x.rstrip(' ').strip(' ') for x in seq_discard_tags.split(',')]
#/
# convert percentage "seq_len_deviation" into fraction
seq_len_deviation = seq_len_deviation/100
#/
##/

## Validate
# check if download workdir exist from previous script
download_workdir = workdir+'/'+'virus_download_wd'
if not os.path.exists(download_workdir):
    print('Was unable to locate download-workdir for viral sequences. Make sure the virus metadata module was run')
    print('Terminating')
    sys.exit()
#/
# check if virus family fasta file exist from previous script
family_fasta_file_path = workdir+'/'+'virus_download_wd'+'/'+'download.fa'
if not os.path.exists(family_fasta_file_path):
    print('Was unable to locate fasta file for viral sequences. Make sure the virus metadata module was run')
    print('Terminating')
    sys.exit()
#/
# check if metadata-file exist from previous script
metadata_file_path = workdir+'/'+'virus_download_wd'+'/'+'metadata_base.tsv'
if not os.path.exists(metadata_file_path):
    print('Was unable to locate metadata file for viral sequences. Make sure the virus metadata module was run')
    print('Terminating')
    sys.exit()
#/
# check so user input any tax to select
if not (tax_target_genus or tax_target_species or tax_target_serotype or tax_target_custom):
    print('No taxonomy selected. Please select a taxonomy using parameters under "Taxonomy Selection". See --help for info')
    print('To overview downloaded taxonomy metadata, see TSV-files at the specified workdir')
    print('Terminating')
    sys.exit()
#/
# Rough check so that "tax_target_custom" has some info to parse
if tax_target_custom and not str(tax_target_custom).count(':') > 0:
    print('Flag --tax_target_custom specified, however, the format could not be interpreted')
    print('Please make sure it is formatted as "key:target" with multiple entries separated by a comma and try again!')
    print('Terminating')
    sys.exit()
#/
##/
###/

### Import data from metadata-module
print('Importing metadata...',flush=True)
## Import metadata file ("base"-info [info from fasta headers])
taxid_metadata_base = {}
with open(metadata_file_path,'r') as f:
    # HEADER: ['taxid','name','num','seq_min','seq_max','seq_med','seq_mean','description']
    for enum,line in enumerate(f):
        if enum == 0: continue # skip header-line
        line = line.strip('\n')
        taxid,name,num,seq_min,seq_max,seq_med,seq_mean,description = line.split('\t')
        num,seq_min,seq_max,seq_med,seq_mean = map(int,[num,seq_min,seq_max,seq_med,seq_mean])
        
        # save
        if taxid in taxid_metadata_base:
            print('WARNING: already imported data for taxid '+taxid)
            print('Terminating!')
            sys.exit()
        taxid_metadata_base[taxid] = {'name':name,'num':num,'seq_min':seq_min,'seq_max':seq_max,'seq_med':seq_med,'seq_mean':seq_mean,'description':description}
        #/
##/
## Import NCBI metadata file (taxname data per taxid)
ncbi_metadata = {}
if not no_ncbi:
    ncbi_metadata_file = download_workdir+'/'+'metadata_ncbi.tsv'
    # confirm that pickle exists
    if not os.path.exists(ncbi_metadata_file):
        print('FATAL: Could not locate NCBI metadata-file. Did you run the metadata-step with option to use NCBI metadata?')
        print(ncbi_metadata_file)
        print('Terminating!')
        sys.exit()
    #/
    # HEADER : ['taxid', 'name', 'num_with_tag', 'superkingdom', 'clade', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'serotype', 'no rank', 'undefined_strain', 'superkingdom_taxid', 'clade_taxid', 'kingdom_taxid', 'phylum_taxid', 'class_taxid', 'order_taxid', 'family_taxid', 'genus_taxid', 'species_taxid', 'serotype_taxid', 'no rank_taxid', 'undefined_strain_taxid']
    header = None
    with open(ncbi_metadata_file,'r') as f:
        for enum,line in enumerate(f):
            # parse line
            line = line.strip('\n')
            line = line.split('\t')
            #/
            # skip header line
            if enum == 0:
                header = line
                continue
            #/
            #/
            # Get taxid+taxname for current entry
            taxid,taxname = line[:2]
            num_genomes_With_tag = line[2]
            #/
            # parse pairs of taxname + taxid: expect X columns first of taxnames followed by X columns of taxids:
            # family,genus,species,family_taxid,genus_taxid,species_taxid
            taxnames_taxids_columns = line[3:]
            taxnames_taxids_header = header[3:]
            number_of_pairs = int(len(taxnames_taxids_columns)/2)
            taxnames_columns = taxnames_taxids_columns[:number_of_pairs] # get taxnames (i.e. adenoviridae, adenovirus, adenovirus B)
            taxnames_columns_header = taxnames_taxids_header[:number_of_pairs] # get header for taxname columns (i.e. family, genus, species)
            taxids_columns = taxnames_taxids_columns[number_of_pairs:] # get taxids (i.e. 1337, 13377, 133777)
            #/
            ## save
            # check so this entry was not saved before
            if taxid in ncbi_metadata:
                print('WARNING: already imported this taxid. Do not expect multiple entries')
                print(taxid)
                print('Terminating!')
                sys.exit()
            #/
            # compile save
            tmp_save = [] # array of [taxlevelname,taxname,taxid]. do array to keep order of taxlevelnames.
            for enum,taxlevelname in enumerate(taxnames_columns_header):
                tmp_save.append( {'taxlevelname':taxnames_columns_header[enum],
                                  'taxname':taxnames_columns[enum],
                                  'taxid':taxids_columns[enum]}
                                )
            #/
            # save
            ncbi_metadata[taxid] = tmp_save
            #/
##/
## Import NCBI metadata pickle-file (taxnames_data in metadata-module)
taxnames_ncbi_data = {} # taxname -> data {taxid, taxlevelname, bvbrc_ids}
taxids_to_taxnames = {} # taxid -> taxname
if not no_ncbi:
    ncbi_metadata_pickle = download_workdir+'/'+'metadata_ncbi_taxnames_data.pickle'
    # confirm that pickle exists
    if not os.path.exists(ncbi_metadata_pickle):
        print('FATAL: Could not locate NCBI metadata pickle-file. Did you run the metadata-step with option to use NCBI metadata?')
        print(ncbi_metadata_pickle)
        print('Terminating!')
        sys.exit()
    #/
    # import ncbi metadata pickle
    with open(ncbi_metadata_pickle,'rb') as f:
        taxnames_ncbi_data = pickle.load(f)
    #/
    # restructure: taxid -> taxname
    for taxname,data in taxnames_ncbi_data.items():
        taxid = data['taxid']
        taxids_to_taxnames[taxid] = taxname
    #/
##/
###/

### Open fasta, filter sequences based on user taxonomy input
def check_taxlevel_entry(taxid,target_level,match_input):
    match_found = False
    if taxid in ncbi_metadata:
        for taxlevel_entry in ncbi_metadata[taxid]:
            if ( (taxlevel_entry['taxlevelname'] == target_level) and 
                (taxlevel_entry['taxname'].lower() in match_input or taxlevel_entry['taxid'] in match_input)
               ):
                   match_found = True
    return match_found

seqheaders_pass_metadata = []
seqheaders_pass_length = []
seqheaders_selected = []
num_parsed = 0
seqs_selected1 = {} # taxid -> bvbrc_id -> seq len
seqs_headers = {} # taxid -> bcbrc_id -> fasta header
print('Reading fasta-sequences, determining which ones to keep based on input taxonomy')
print('target_genus '+str(tax_target_genus))
print('target_species '+str(tax_target_species))
print('target_serotype '+str(tax_target_serotype))
print('target_custom '+str(tax_target_custom),flush=True)
for entry in SeqIO.parse(family_fasta_file_path,'fasta'):
    ## Get entry info
    # parse description
    description = entry.description
    #/
    # parse bvbrc_id
    bvbrc_id = description.split('| ')[1].split(']')[0] # get BVBRC ID from fasta header. Format:  'description': '5KW1_C   Chain C, DNA/RNA (30-MER).   [unidentified adenovirus | 10535.961]'
    #/
    # get taxid from bvbrc id
    taxid = bvbrc_id.split('.')[0]
    #/
    # get length of sequence
    seq_len = len(entry.seq)
    #/
    ##/
    
    ### "Base" metadata checks
    ## Run checks
    if 0 and 'TODO':
        pass
    ##/
    ## Check if sequence is to be included
    ##/
    ###/
    
    ### NCBI metadata checks
    ## Run checks
    inclusion_checks = []
    # genus
    if tax_target_genus != None:
        if check_taxlevel_entry(taxid,'genus',tax_target_genus):
            inclusion_checks.append('match_genus')
    #/
    # species
    if tax_target_species != None:
        if check_taxlevel_entry(taxid,'species',tax_target_species):
            inclusion_checks.append('match_species')
    #/
    # serotype
    if tax_target_serotype != None:
        if check_taxlevel_entry(taxid,'serotype',tax_target_serotype):
            inclusion_checks.append('match_serotype')
    #/
    # custom
    if tax_target_custom != None:
        for match_key_and_val in tax_target_custom: # expected format: array of [key1:val1,key2:val2]
            match_key,match_val = match_key_and_val.split(':')
            if check_taxlevel_entry(taxid,match_key,[match_val]):
                inclusion_checks.append('match_custom')
    #/
    ##/
    ## Check if sequence is to be included
    failed_target_checks = 0
    if inclusion_checks:
        if tax_target_genus != None and not 'match_genus' in inclusion_checks:      failed_target_checks += 1
        if tax_target_species != None and not 'match_species' in inclusion_checks:      failed_target_checks += 1
        if tax_target_serotype != None and not 'match_serotype' in inclusion_checks:      failed_target_checks += 1
        if tax_target_custom != None and not 'match_custom' in inclusion_checks:      failed_target_checks += 1
        
        if failed_target_checks == 0:
            # save headerline
            seqheaders_pass_metadata.append(description)
            #/
            # save seq len
            if not taxid in seqs_selected1:              seqs_selected1[taxid] = {}
            if bvbrc_id in seqs_selected1[taxid]: print('WARNING: already had a length saved at this BVBRC ID')
            seqs_selected1[taxid][bvbrc_id] = seq_len
            #/
            # save description at taxid->bvbrc_id
            if not taxid in seqs_headers:           seqs_headers[taxid] = {}
            seqs_headers[taxid][bvbrc_id] = description
            #/
    ##/
    ###/
    num_parsed += 1

print('After filtering input matches there are '+str(len(seqheaders_pass_metadata))+' sequences',flush=True)
###/

### Terminate script if there are 0 sequences passed
if not len(seqheaders_pass_metadata) > 0:
    print('No sequences were parsed. Please reconsider your selection options')
    print('Terminating!')
    sys.exit()
###/

### Determine length-selection
print('Determining range of sequence lengths to parse',flush=True)
target_len_range = None
if seq_len_range: # user already defined a range of lengths to parse
    target_len_range = seq_len_range
else:
    # determine which sequences to use for length-range estimate
    INFO_taxids_without_metadata = set()
    INFO_taxids_with_metadata = set()
    lens_statistica = [] # store median/mean lengths
    lens_tot = [] # store min/max reported lengths
    for taxid in seqs_selected1:
        # check if this taxid had no ncbi "base" metadata (not expected)
        if not taxid in taxid_metadata_base:
            if 0 and 'verbose':
                print('INFO: Did not find taxid in NCBI "base" metadata. This is not expected. Taxid: '+str(taxid))
                print('Will skip metadata for this taxid and proceed like this did not happen...')
            INFO_taxids_without_metadata.add(taxid)
            continue
        #/
        # get "base" metadata
        taxid_metadata = taxid_metadata_base[taxid]
        INFO_taxids_with_metadata.add(taxid)
        #/
        # Add length "statistica"
        if seq_len_method == 'median_of_medians':
            lens_statistica.append(taxid_metadata['seq_med'])
        elif seq_len_method == 'mean_of_means':
            lens_statistica.append(taxid_metadata['seq_mean'])
        #/
        # Add length min/max
        if seq_len_select_all:
            lens_tot.append(taxid_metadata['seq_min'])
            lens_tot.append(taxid_metadata['seq_max'])
        #/
    #/
    # print user info
    if INFO_taxids_without_metadata:
        print('INFO: Did not find taxid in NCBI "base" metadata for '+str(len(INFO_taxids_without_metadata))+' taxids. Found metadata for '+str(len(INFO_taxids_with_metadata))+' taxids')
        print('If the number of taxids without NCBI "base" metadata is very (>90%) high then this indicates an optimal length-representative for all taxids is not available. You may need to specify a custom range for length using --length_range')
    #/
    # determine length-range
    if seq_len_select_all:
        target_len_range = [min(lens_tot),max(lens_tot)]
    else:
        # determine which value to use as "midpoint" in the distribution of lengths
        len_statistica2 = None
        if seq_len_method == 'median_of_medians':
            len_statistica2 = median(lens_statistica)
        elif seq_len_method == 'mean_of_means':
            len_statistica2 = mean(lens_statistica)
        #/
        # compile the range with deviations around the "midpoint"
        target_len_range = [ int(len_statistica2*(1-seq_len_deviation)) , int(len_statistica2*(1+seq_len_deviation)) ]
    #/
print('Range of lengths determined at: '+str(target_len_range[0])+'-'+str(target_len_range[-1]),flush=True)
###/

### Determine which sequences to parse
seqs_selected2 = {} # taxid -> bvbrc_id -> seq_len
num_passed_len = 0
num_passed_tag = 0
num_passed = 0
for taxid,bvbrc_ids_lens in seqs_selected1.items():
    for bvbrc_id,seq_len in bvbrc_ids_lens.items():
        # check if pass len
        pass_len = False
        if seq_len >= target_len_range[0] and seq_len <= target_len_range[-1]:
            pass_len = True
            num_passed_len += 1
        #/
        # check discard tag is supplied and sequence pass it
        pass_discard_tag = True # assume it is true, prove the entry wrong
        if seq_discard_tags != None:
            for discard_tag in seq_discard_tags:
                if seqs_headers[taxid][bvbrc_id].find(discard_tag) != -1:
                    pass_discard_tag = False
        if pass_discard_tag and seq_discard_tags != None:
            num_passed_tag += 1
        #/
        
        # check if save
        if pass_len and pass_discard_tag:
            if not taxid in seqs_selected2:     seqs_selected2[taxid] = {}
            seqs_selected2[taxid][bvbrc_id] = seq_len
            num_passed += 1
        #/

print('After filtering sequence lengths there are '+str(num_passed)+' sequences (by length: '+str(num_passed_len)+', by discard-tag [if supplied]: '+str(num_passed_tag)+')',flush=True)
###/

### Check if user wants some info of selection
## Check if user wants to print the fasta headers
if perform_fasta_header_print:
    for taxid in seqs_selected2:
        for bvbrc_id in seqs_selected2[taxid]:
            print(seqs_headers[taxid][bvbrc_id])
##/
## Check if user wants to "glance" the selection results. Then quit here
if perform_glance:
    print('Argument --glance specified, terminating now!')
    sys.exit()
###/

### Parse sequences
# init genomes-folder
genomes_dir = workdir+'/'+'genomes'
if os.path.exists(genomes_dir):
    print('Previous directory for genomes found. Will remove it before populating it with new sequences!')
    shutil.rmtree(genomes_dir)
os.makedirs(genomes_dir)
#/
# Traverse fasta-file and output sequences that passed
print('Writing genome sequences',flush=True)
num_genome_sequences_written = 0
for entry in SeqIO.parse(family_fasta_file_path,'fasta'):
    ## Get entry info
    # parse description
    description = entry.description
    #/
    # parse bvbrc_id
    bvbrc_id = description.split('| ')[1].split(']')[0] # get BVBRC ID from fasta header. Format:  'description': '5KW1_C   Chain C, DNA/RNA (30-MER).   [unidentified adenovirus | 10535.961]'
    #/
    # get taxid from bvbrc id
    taxid = bvbrc_id.split('.')[0]
    #/
    # determine file name (first entry after splitting header. I.e. SKW1_C in above example)
    fasta_name = description.split()[0]
    #/
    ##/
    
    ## Check if sequence pass and write it
    if taxid in seqs_selected2 and bvbrc_id in seqs_selected2[taxid]:
        # determine file name
        fasta_file_name = genomes_dir+'/'+fasta_name+'.fasta'
        #/
        # warn user if this fasta file already existed
        if os.path.exists(fasta_file_name):
            print('WARNING: This fasta file was already written. Sequences are expected to have an unique ID in their headers')
            print(fasta_file_name)
            
            if not fasta_duplicates_do_not_terminate:
                print('Apply flag --ignore_duplicates to ignore this warning and proceed, will terminate now.')
                print('Terminating!')
                sys.exit()
        #/
        # write sequence to file
        with open(fasta_file_name,'w') as nf:
            nf.write('>'+description+'\n'+str(entry.seq)+'\n')
            num_genome_sequences_written += 1
    ##/
#/
print('Number of genomes written: '+str(num_genome_sequences_written))
print('Done!',flush=True)
###/

### Check if user wants to generate an outgroup
if not do_not_generate_outgroup:
    print('Begin outgroup determination',flush=True)
    
    ## Remove previous outgroup, if it exists
    # outgroup dir with genome file
    if os.path.exists(workdir+'/'+'outgroup'):
        print('A previously generated directory for outgroup found, will remove before proceeding...')
        shutil.rmtree(workdir+'/'+'outgroup')
    #/
    # outgroup file of outgroup header
    if os.path.exists(workdir+'/'+'outgroup_accession.txt'):
        print('A previously generated file for outgroup found, will remove before proceeding...')
        os.remove(workdir+'/'+'outgroup_accession.txt')
    #/
    ##/
    
    ## Init WD to determine outgroup
    outgroup_wd = workdir+'/'+'virus_outgroup_wd'
    if os.path.exists(outgroup_wd):
        print('Previous workdir for outgroup identified, removing it before creating a new one')
        shutil.rmtree(outgroup_wd)
    os.makedirs(outgroup_wd)
    ##/

    ## Get taxids which are outgroup candidates
    outgroup_candidate_taxids = set()
    for taxid,base_data in taxid_metadata_base.items():
        # skip this taxid if it is part of selected output
        if taxid in seqs_selected2: continue
        #/
        # select only other taxids where at least X number of other datasets are available
        if base_data['num'] < outgroup_candidates_taxid_min_genomes: continue
        #/
        # save
        outgroup_candidate_taxids.add(taxid)
        #/
    
    if not outgroup_candidate_taxids:
        print('Did not find any datasets to use as outgroups. Make sure your selection does not include the whole virus family. If this issue still persists then lower the "--outgroup_candidates_taxid_min_genomes" or set it to 0')
        print('Terminating!')
        sys.exit()
    ##/
    ## Parse sequences and output to temporary genomes directory
    print('Adding sequences to temporary genomes directory',flush=True)
    outgroup_genomes_dir = outgroup_wd+'/'+'genomes'
    os.makedirs(outgroup_genomes_dir)
    
    max_per_taxid = 3
    num_written_per_taxid = {}
    for entry in SeqIO.parse(family_fasta_file_path,'fasta'):
        ## Get entry info
        # parse description
        description = entry.description
        #/
        # parse bvbrc_id
        bvbrc_id = description.split('| ')[1].split(']')[0] # get BVBRC ID from fasta header. Format:  'description': '5KW1_C   Chain C, DNA/RNA (30-MER).   [unidentified adenovirus | 10535.961]'
        #/
        # get taxid from bvbrc id
        taxid = bvbrc_id.split('.')[0]
        #/
        # get length of sequence
        seq_len = len(entry.seq)
        #/
        ##/
        ## Run checks
        # skip if we already wrote X number of sequences for this taxid
        if taxid in num_written_per_taxid and num_written_per_taxid[taxid] > max_per_taxid: continue
        #/
        # skip if not part of selected outout or outgroup candidates
        if not (taxid in seqs_selected2 or taxid in outgroup_candidate_taxids): continue
        #/
        # skip if not within certain range of "midpoint" sequence length metric
        if not taxid in taxid_metadata_base: continue
        seqlen_metric = taxid_metadata_base[taxid]['seq_med']
        if not (seq_len >= seqlen_metric*0.9 and seq_len <= seqlen_metric*1.1): continue
        #/
        ##/
        ## Write sequence
        prefix = None
        if taxid in seqs_selected2:                 prefix = 'S' # add "S" for "selected" identifier for sequences that exist in selection
        if taxid in outgroup_candidate_taxids:      prefix = 'O' # add "O" for "outgroup" identifier for sequences that are outgroup candidates
        
        fasta_name = description.split()[0]
        fasta_file_name = outgroup_genomes_dir+'/'+prefix+'_'+fasta_name+'.fasta'
        with open(fasta_file_name,'w') as nf:
            nf.write('>'+description+'\n'+str(entry.seq)+'\n')
        
        if not taxid in num_written_per_taxid:      num_written_per_taxid[taxid] = 0
        num_written_per_taxid[taxid] += 1
        ##/
    ##/
    ## Execute mashtree
    print('Executing mashtree',flush=True)
    mashtree_cmd = ['mashtree',
                    '--genomesize',int(mean(target_len_range)),
                    '--mindepth','0',
                    #'--sketch-size',1000, ## the default sketch-size [10 kbp] empirically make "best" results (distance matrixs has fewer values at 1)
                    '--outmatrix',outgroup_wd+'/'+'distance_matrix.mashtree.tsv',
                    outgroup_genomes_dir+'/'+'*.fasta','>',outgroup_wd+'/'+'mashtree.dnd']
    subprocess.call(' '.join(map(str,mashtree_cmd)),shell=True)
    if 0 and 'VERBOSE':
        print('Mashtree CMD:\n'+' '.join(map(str,mashtree_cmd)))
    ##/
    ## Validate that mashtree executed with output
    if not os.path.exists(outgroup_wd+'/'+'distance_matrix.mashtree.tsv'):
        print('Could not locate Mashtree output. Make sure that it is installed and can be executed')
        print('Terminating!')
        sys.exit()
    ##/
    ## Parse distance-matrix
    print('Parsing mashtree distance-matrix',flush=True)
    # parse file
    header = None
    row_data = []
    with open(outgroup_wd+'/'+'distance_matrix.mashtree.tsv','r') as f:
        for enum,line in enumerate(f):
            # parse line
            line = line.strip('\n')
            line = line.split('\t')
            #/
            # check if header line, then save
            if enum == 0:
                header = line
                continue
            #/
            # parse row
            row_data.append(line)
            #/
    #/
    # compute all-vs-all comparisons
    all_vs_all_comps = {}
    for idx,seq_name in enumerate(header):
        if idx == 0: continue # skip column for "dataset row-identifier"
        for row in row_data:
            oSeq_name = row[0]
            # check if comp is vs. self
            if seq_name == oSeq_name: continue
            #/
            # get val
            val = float(row[idx])
            #/
            # determine "comparison ID"
            comp_ID = '||'.join(sorted([seq_name,oSeq_name]))
            #/
            # save
            all_vs_all_comps[comp_ID] = val
            #/
    #/
    # compute intra-S, intra-O, and S-vs-O values
    groups_comps = {}
    for comp_ID,val in all_vs_all_comps.items():
        # determine what comp this is
        seq1,seq2 = comp_ID.split('||')
        SK = None

        if seq1[0] == 'S' and seq2[0] == 'S':       SK = 'S-vs-S' # check if S-vs-S
        elif seq1[0] == 'O' and seq2[0] == 'O':       SK = 'O-vs-O' # check if O-vs-O
        elif (seq1[0] == 'S' and seq2[0] == 'O') or (seq1[0] == 'O' and seq2[0] == 'S'):       SK = 'S-vs-O' # check if S-vs-O
        #/
        # save val
        if not SK in groups_comps:      groups_comps[SK] = []
        groups_comps[SK].append(val)
        #/
    #/
    # compute the median distance within groups
    groups_comps_stats = {}
    for group,values in groups_comps.items():
        groups_comps_stats[group] = {'median':median(values),'mean':mean(values),'stdev':stdev(values)}
    #/
    # print some info
    print('Distances between groups (selected vs. selected, outgroup_candidate vs. outgroup_candidate, selected vs. outgroup_candidate')
    print('group\tmedian\tmean  \tstdev:')
    for group,stats in sorted(groups_comps_stats.items(),key=lambda x: x[1]['median']):
        printArr = [group,round(stats['median'],4),round(stats['mean'],4),round(stats['stdev'],4)]
        print('\t'.join(map(str,printArr)))
    print('If the value for S-vs-S is NOT the lowest then there is a risk that the outgroup is poorly assigned',flush=True)
    
    if sorted(groups_comps_stats.items(),key=lambda x: x[1]['median'])[0][0] != 'S-vs-S':
        print('WARNING: S-vs-S did not have the lowest distance. You might want to manually specify an outgroup')
        print('Proceeding like this did not happen',flush=True)
    #/
    ##/
    ## Select outgroup (take the outgroup candidate ["O"] with the lowest maximum-value to selected ["S"]) 1 = "very distant" 0 = "identical"
    # get all O_seq values per S_seq
    S_vs_O_values = {} # S_seqX -> [O_seq values]
    O_vs_S_values = {} # O_seqX -> [S_seq values]
    for idx,seq1 in enumerate(header):
        if idx == 0: continue # skip column for "dataset row-identifier"
        # skip if seq1 is not a selected "S" sequence
        if not seq1[0] == 'S': continue
        #/
        for row in row_data:
            seq2 = row[0]
            # skip if seq2 is not an outgroup candidate "O" sequence
            if not seq2[0] == 'O': continue
            #/
            # get val
            val = float(row[idx])
            #/
            # save
            if not seq1 in S_vs_O_values:       S_vs_O_values[seq1] = {}
            S_vs_O_values[seq1][seq2] = val
            #/
            # save2
            if not seq2 in O_vs_S_values:       O_vs_S_values[seq2] = {}
            O_vs_S_values[seq2][seq1] = val
            #/
    #/
    # get maximum value per S-seq
    Oseq_to_Sseq_minmaxVals = {}
    for Oseq,Sseq_vals in O_vs_S_values.items():
        Sseq_vals_sorted = sorted(Sseq_vals.items(),key=lambda x: x[1])
        minval = Sseq_vals_sorted[0]
        maxval = Sseq_vals_sorted[-1]
        Oseq_to_Sseq_minmaxVals[Oseq] = {'min':minval,'max':maxval}
    #/
    # sort by best candidate
    Oseq_best_candidates = sorted(Oseq_to_Sseq_minmaxVals.items(),key=lambda x: (-x[1]['max'][1],x[1]['min'][1]))
    #/
    # select the candidate which exist at "median" S-vs-O group: Rationale is that this should give enough margin to "always" be outside the S-cluster in a phylogenetic tree.
    # If we take the BEST candidate based on the distance matrix, in theory we should get the O-dataset that is closest to S-cluster.
    # However, empirically the S-cluster not always perfectly separated from O-sequences.
    outgroup_candidate_selected = None
    for candidate,distance_data in Oseq_best_candidates:
        # determine a value to test against to select the optimal outgroup dataset
        test_against_value = (groups_comps_stats['S-vs-O']['mean']-groups_comps_stats['S-vs-O']['stdev']) # default: take first value that is greater than the average distance minus one standard deviation
        if test_against_value < groups_comps_stats['S-vs-O']['median']:     test_against_value = groups_comps_stats['S-vs-O']['median'] # however, if that value is very small (indicative of heterogeneous O-datasets), use the median instead.
        #/
        # run test
        if distance_data['min'][1] >= test_against_value:
            outgroup_candidate_selected = [candidate,distance_data]
            break # break on first
        #/
    #/
    # check if any outgroup was determined
    if outgroup_candidate_selected == None:
        print('WARNING: Unable to assign an outgroup.')
        print('Terminating!')
        sys.exit()
    #/
    # print info
    print('Determined outgroup sample:')
    print(outgroup_candidate_selected[0] + ' at distance (0-1): min_highest='+str(round(outgroup_candidate_selected[1]['min'][1],4))+' max_highest='+str(round(outgroup_candidate_selected[1]['max'][1],4))+' distance')
    #/
    ##/
    ## Save outgroup to outer
    print('Writing outgroup',flush=True)
    # get identifier of selected candidate
    outgroup_identifier = outgroup_candidate_selected[0]
    if outgroup_identifier[:2] == 'O_':      outgroup_identifier = outgroup_identifier[2:]
    #/
    # Write sequence identifier into "<wd>/outgroup.txt"
    with open(workdir+'/'+'outgroup_accession.txt','w') as nf:
        nf.write(outgroup_identifier+'\n')
    #/
    # Put sequence into folder "<wd>/outgroup/<file.fasta>"
    outgroup_genome_dir = workdir+'/'+'outgroup'
    os.makedirs(outgroup_genome_dir)
    
    with open(outgroup_genome_dir+'/'+outgroup_identifier+'.fasta','w') as nf:
        for entry in SeqIO.parse(family_fasta_file_path,'fasta'):
            # parse description
            description = entry.description
            #/
            # if description has the selected outgroup id, then write
            if description.split()[0] == outgroup_identifier:
                nf.write('>'+description+'\n'+str(entry.seq)+'\n')
            #/
    #/
    ## check if clear temporary folder
    if not keep_files:
        print('Clearing temporary files')
        shutil.rmtree(outgroup_wd)
    ##/
    
    print('Number of genomes written: '+str(num_genome_sequences_written)+', outgroup: 1')
    print('Done!',flush=True)
    ##/
###/