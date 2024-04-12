#!/usr/bin/env python3

import os
import sys
import argparse
import pickle
import shutil
from Bio import SeqIO
from statistics import mean,median
from time import sleep
import requests

## Static variables
# bvbrc ftp details
bvbrc_ftp = 'ftp.bvbrc.org'
bvbrc_ftp_user = 'anonymous'
bvbrc_ftp_password = ''
bvbrc_ftp_file_path = 'viruses'
#/
# used for NCBI metadata: the order in which to present taxonomic names (the last entry should always be "undefined_strain" which is assigned to taxids that do not match at any of NCBI levels)
taxnames_to_parse_ordered = ['superkingdom','clade','kingdom','phylum','class','order','family','genus','species','serotype','no rank']
undefined_strain_identifier = 'undefined_strain'
taxnames_to_parse_ordered.append(undefined_strain_identifier)
#/
##/

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = 
'''
Retrieve data from BV-BRC for selected taxonomy. Virus-equivalent for "metadata" module
''')

parser.add_argument('-wd','--workdir',required=True,help='Path to working directory. This folder will be created if not already present')

parser.add_argument('-tf','--target_family',required=True,help='Target family (Example: adenoviridae)')

parser.add_argument('-f','--filter',required=False,default='complete genome',help='Filter for fasta headers to use when determining sequences to use for length-estimate (default: use sequences with headers containing "complete genome" to determine length)')

#parser.add_argument('-ncbi','--ncbi_metadata',required=False,action='store_true',default=True,help='[REQUIRED AS OF NOW FOR DOWNSTREAM] If specified, will download additional metadata (e.g taxonomic names) from NCBI using Entrez (default: do not download NCBI metadata)')

#parser.add_argument('--limit',type=int,default=None,help='Limits the analysis to the specified number of genomes (e.g. for test/debug purposes)')
#/
# parse input
if 1 and 'run':
    args = parser.parse_args()
else:
    print('IDE MODE')
    args = parser.parse_args(['--workdir','IDE_virus',
                              '--target_family','Adenoviridae',
                              '--ncbi_metadata'])

target_family = args.target_family

workdir = args.workdir

filter_require_fasta_match = args.filter

parse_metadata_from_ncbi = True
#/
# validate input
if not target_family:
    print('Must supply a family as input using --target_family (e.g., adenoviridae)')
    sys.exit()
    
#/
# Final formatting
if target_family:           target_family = target_family.lower()
#/
###/

### Download family fasta
# Init download workdir
download_workdir = workdir+'/'+'virus_download_wd'
if not os.path.exists(download_workdir):        os.makedirs(download_workdir)
#/
## Download family fasta from bvbrc
target_family_capitalized = target_family[0].upper()+target_family[1:] # need to capitalize the first letter to match files at bvbrc ftp (e.g. Adenoviridae)
download_file_path = download_workdir+'/'+'download.fa'
if not (os.path.exists(download_file_path) and os.path.getsize(download_file_path) > 0):
    # download using python ftplib
    print('Need to download file from bvbrc, executing download now (using anonymous/anonymous at ftp login)',flush=True)
    from ftplib import FTP
    with FTP(bvbrc_ftp) as ftp_handle:
        ftp_handle.login(user='anonymous', passwd='')
        
        # define the target file path at ftp (e.g. <ftp_base_url>/viruses/Virus.fna)
        ftp_target_file_path = bvbrc_ftp_file_path+'/'+target_family_capitalized+'.fna'
        #
        # Get size of file at remote ftp (will also test that the remove file exists, i.e. the user input a valid name)
        ftp_handle.sendcmd('TYPE I') # Set binary mode in ftp_handle to allow for "size"-call (ascii mode does not allow that)
        try:
            remote_file_size = ftp_handle.size(ftp_target_file_path)
        except:
            print('WARNING: Could not get the size of remote file for input '+target_family_capitalized)
            print('Please make sure that this family exists at BV-BRC!')
            print('NOTE: some virus families exist in the collection "viruses" at BV-BRC, for example Arenaviridae (as of writing this)')
            print('You can download these genomes using --target_family viruses')
            print('Your target family (or lower level lineages) are parsed in the vgenome module')
            print('Terminating!')
            sys.exit()
        #/
        # Try to download the file from the ftp
        ftp_handle.sendcmd('TYPE A') # Set ASCII-mode to allow non-binary download for convenient parsing
        
        def write_lines(inp_data): # need a function to add newlines to each line parsed from ftp. apparently retrlines in ftp does not send newline characters, but it does send each line as one batch.
            f.write(inp_data+'\n')
            
        with open(download_file_path, 'w') as f:
            ftp_handle.retrlines('RETR ' + ftp_target_file_path, write_lines)
        #/
        # Get size of downloaded file
        download_file_size = 0
        if os.path.exists(download_file_path):
            download_file_size = os.path.getsize(download_file_path)
        #/
        # Check if sizes match
        print('Verifying download...',flush=True)
        if remote_file_size == download_file_size:
            print('Local and remote files have the same size [OK]',flush=True)
        elif abs(remote_file_size-download_file_size) < 1000:
            print('Local and remote files diff less than 1000 bytes [OK]',flush=True)
        else:
            print('WARNING: local and remote files difference is greater than expected:')
            print('Remote: ' +str(remote_file_size))
            print('Local: ' +str(download_file_size))
            print('Terminating!')
            sys.exit()
        #/
        
    print('Download finished',flush=True)
    #/
else:
    print('File already downloaded from bvbrc, proceeding',flush=True)
##/
###/

### Compile metadata from downloaded fasta (+optionally from other sources)
def get_taxon_data_from_entrez(tax_ids_list,num_ids_per_query=100):
    def send_entrez_query(taxon_ids):
        # send request to entrez and get response (manual: https://www.ncbi.nlm.nih.gov/books/NBK25499/; see table 1)
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        entrez_params = {
                        'db': 'taxonomy',
                        'id': taxon_ids,
                        'retmode': 'xml'
                  }
        response = requests.get(base_url + 'efetch.fcgi', params=entrez_params)
        response_content = response.text
        response_content2 = response_content.split('\n') # return as lists of lines
        return response_content2
        #/
    ## Convert input to list if it is not already
    if not type(tax_ids_list) == list:
        tax_ids_list = list(tax_ids_list)
    ##/
    
    ## Split input tax_id_list into lists of size X and submit them as queries to entrez. Accumulate the output of batches to a "non-batched"/"single-input" query.
    tax_ids_sublists = [tax_ids_list[i:i+num_ids_per_query] for i in range(0, len(tax_ids_list), num_ids_per_query)] # make sublists of size X
    entrez_response = []
    for enum,sublist in enumerate(tax_ids_sublists):
        print('Submitting sublist '+str(enum)+' to Entrez with '+str(len(sublist))+' taxonomy IDs...',flush=True)
        tmp_entrez_response = send_entrez_query(sublist)
        entrez_response += tmp_entrez_response
        sleep(0.5) # Make a short sleep to not be banned from entrez. It is capped at 3 request per second.
        print('... done!',flush=True)
    ##/
    
    ## chomp response, grouped by input taxid. Input: lines from response, output: groups as strings of lines
    # the structure of return is:
    # <Taxon> # taxid1
    #      data
    # </Taxon>
    # <Taxon> # taxid2
    #      data
    # </Taxon>
    response_grouped = []
    for line in entrez_response:
        # skip empty lines
        if not line: continue
        #/
        # check if init new group (intendation is 0 and <Taxon> identifier is present)
        if line[0] != ' ' and line.find('<Taxon>') != -1:
            response_grouped.append('')
        #/
        # check if row has intendation, then add the content to current group
        elif line[0] == ' ':
            response_grouped[-1] += line
        #/
    ##/
    
    ## for each parsed taxid, get data
    taxids_data = {} # taxid -> entrez data
    for chunk_raw in response_grouped:
        # Get content before info about parents etc in "tree of life"
        chunk = chunk_raw.split('<LineageEx>')[0]
        #/
        # Get taxid
        taxid = chunk.split('<TaxId>')[1].split('</TaxId>')[0]
        #/
        # Skip taxid if not part of input. Dont know why we get back something that we didnt request
        if not taxid in tax_ids_list:
            print('Got an unexpected taxid in return from Entrez: '+str(taxid)+' Ignoring it...',flush=True)
            continue
        #/
        # Get scientific name
        scientific_name = chunk.split('<ScientificName>')[1].split('</ScientificName>')[0]
        #/
        # Get rank
        rank = chunk.split('<Rank>')[1].split('</Rank>')[0]
        #/
        # IDE: see if there is genome length information
        if chunk.find('genome') != -1:
            print('IDE: genome-information found: ')
            print(chunk)
        #/
        # Get info from organisms upstream in tree of life
        tree_of_life_data = chunk_raw.split('<LineageEx>')[1].split('</LineageEx>')[0]
        level_data = {} # taxonomic level -> data
        level_data[rank] = {'taxid':taxid,'name':scientific_name,'level':rank} # If current input taxid is not far down in "tree of life" then save its info. E.g., taxid at family-level only contain info in its tree upstreams, not downstreams.
        for chunk2_raw in tree_of_life_data.split('</Taxon>'):
            # check if out of boundaries. expect Taxon-tag to exist
            if chunk2_raw.find('<Taxon>') == -1: continue
            #/
            chunk2 = chunk2_raw.split('<Taxon>')[1] # make sure there is no residue content prior to <Taxon> entry
            chunk_taxid = chunk2.split('<TaxId>')[1].split('</TaxId>')[0]
            chunk_scientific_name = chunk2.split('<ScientificName>')[1].split('</ScientificName>')[0]
            chunk_tree_level = chunk2.split('<Rank>')[1].split('</Rank>')[0]
            
            # bugcheck
            if chunk_tree_level in level_data and not chunk_tree_level == 'no rank': # for some reason, "no rank" is stated in the tree but no other taxonomic level. So for this rankname we allow duplicate entries
                print('WARNING: had multiple entries for tree level: '+chunk_tree_level+' taxid: '+taxid)
                print('... Pretending like this did not happen and moving on the with the latest entry')
            #/
            
            # print a warning if we do not have this taxonomic level defined in "taxnames_to_parse_ordered" (global variable)
            if not chunk_tree_level in taxnames_to_parse_ordered:
                print('WARNING: found a taxonomic level of which the order was not determined: '+chunk_tree_level)
                print('Current order: '+', '.join(taxnames_to_parse_ordered))
                print('... Pretending like this did not happen and moving on')
            #/
            # save
            level_data[chunk_tree_level] = {'taxid':chunk_taxid,'name':chunk_scientific_name,'level':chunk_tree_level}
            #/
        #/
        # add "None"-entries for taxonomic names that were not present at current taxid
        for taxname in taxnames_to_parse_ordered:
            if not taxname in level_data:
                level_data[taxname] = {'taxid':None,'name':None,'level':None}
        #/
        # check if current taxid was assigned at any taxonomic level. if not, then save it at the "undefined_strain_identifier"
        taxid_found = False
        for taxlevelname,data in level_data.items():
            if data['taxid'] != None and data['taxid'] == taxid:
                taxid_found = True
                break # break on first occ
        if not taxid_found:
            level_data[undefined_strain_identifier] = {'taxid':taxid,'name':scientific_name,'level':undefined_strain_identifier}
        #/
        # save to outer
        taxids_data[taxid] = {'taxid':taxid,'name':scientific_name,'taxdata':level_data}
        #/
    ##/
    ## check if any of input taxids did not get a return
    taxids_missing_data = set(tax_ids_list).difference(set(taxids_data))
    if taxids_missing_data:
        print('WARNING: Failed to obtain Entrez-data for '+str(len(taxids_missing_data))+' taxonomy IDs:')
        print(','.join(taxids_missing_data))
        print('\nIf this number is high (>100 or >10% of input taxids) then it indicates something went wrong in the request to NCBI and this script should be re-run')
        print('Adding empty entries for these taxon IDs and moving on...',flush=True)
        for taxid in taxids_missing_data:
            tmp_taxdata = {}
            for taxname in taxnames_to_parse_ordered: # taxnames_to_parse_ordered is a global variable
                tmp_taxdata[taxname] = {'taxid':None,'name':None,'level':None}
            taxids_data[taxid] = {'taxid':taxid,'name':None,'taxdata':tmp_taxdata}
    ##/
    return taxids_data,taxids_missing_data

## REVISE DOCSTRINGS. Filter sequences and output them 1 by 1 (BVBRC-ids are strucutred as <ncbi_tax_id>.<enum>)
# Get length of all sequences, determine which ones to keep
seq_lens = {}
taxon_ids_all = set() # save all taxids, even those that do not pass "tag"-filter (e.g. "complete genome")
taxon_ids_all_bvbrc_ids = {} # taxid -> bvbrc_ids
taxon_ids = {} # taxid->count
taxon_ids_randomDescription = {} # taxid -> fasta header
taxon_ids_lens = {} #taxid -> [lens_of_filtered_genomes] # save only the length of genomes with user-defined tag (e.g., "complete genome")
for entry in SeqIO.parse(download_file_path,'fasta'):
    # parse bvbrc_id
    bvbrc_id = entry.description.split('| ')[1].split(']')[0] # get BVBRC ID from fasta header. Format:  'description': '5KW1_C   Chain C, DNA/RNA (30-MER).   [unidentified adenovirus | 10535.961]'
    #/
    # get taxid from bvbrc id
    taxid = bvbrc_id.split('.')[0]
    #/
    # save taxid in "all" collection of names
    taxon_ids_all.add(taxid)
    #/
    # save bvbrc id at taxid
    if not taxid in taxon_ids_all_bvbrc_ids:    taxon_ids_all_bvbrc_ids[taxid] = set()
    taxon_ids_all_bvbrc_ids[taxid].add(bvbrc_id)
    #/
    # check if this sequence has selection indication (i.e. "complete genome")
    if entry.description.find(filter_require_fasta_match) == -1:
        continue
    #/
    # get length of sequence
    seq_len = len(entry.seq)
    #/
    # bugcheck: we do not expect any bvbrc id to exist multiple times in fasta
    if bvbrc_id in seq_lens:
        print('FATAL: BV-BRC id was already processed: '+bvbrc_id)
        print('Terminating!')
        sys.exit()
    #/
    # save len
    seq_lens[bvbrc_id] = seq_len
    #/
    # save number of datasets per taxid
    if not taxid in taxon_ids:          taxon_ids[taxid] = 0
    taxon_ids[taxid] += 1
    #/
    # save fasta header of taxid (indicator-description of a taxid without downloading summary-table or parsing entrez)
    if not taxid in taxon_ids_randomDescription:
        species_description = entry.description.split('[')[1].split(' | ')[0] # get "unidentified adenovirus"-part of description example above
        full_header = entry.description.replace(entry.name,'').lstrip().split('[')[0].strip() # get "Chain C, DNA/RNA (30-MER)."-part of description example above
        taxon_ids_randomDescription[taxid] = {'species':species_description,'full':full_header}
    #/
    # save lengths per taxid
    if not taxid in taxon_ids_lens:     taxon_ids_lens[taxid] = []
    taxon_ids_lens[taxid].append(seq_len)
    #/

print('Identified '+str(len(taxon_ids_all))+' taxon IDs, of which '+str(len(taxon_ids))+' had information according to filtering-tag ("'+filter_require_fasta_match+'") about sequence length in '+str(len(seq_lens))+' sequences:')
print('Max length: '+str(max(seq_lens.values())))
print('Min length: '+str(min(seq_lens.values())))
print('Median length: '+str(int(median(seq_lens.values()))))
print('Mean length: '+str(int(mean(seq_lens.values()))),flush=True)
#/
##/

### Structure metadata for family
## Generate "BASE" metadata (info available from fasta headers)
taxid_metadata_base  = {} # taxid -> data
# init data holder at each taxid
for taxid,num_found in taxon_ids.items():
    taxid_metadata_base[taxid] = {'num':num_found}
#/
# add min / max / median sequence lengths
for taxid,lens in taxon_ids_lens.items():
    taxid_metadata_base[taxid].update( {'min':min(lens),'max':max(lens),'mean':int(mean(lens)),'median':int(median(lens))} )
#/
# add a random description taken from any header of a taxonomy
for taxid,desc in taxon_ids_randomDescription.items():
    taxid_metadata_base[taxid]['desc'] = desc
#/
##/
## Generate "NCBI" metadata (info parsed from ncbi via entrez)
# check if user wants to parse additional data from NCBI entrez
taxid_metadata_ncbi = {}
taxid_metadata_ncbi_namesPerTaxLevel = {} # taxlevel -> nameVariants -> counts_of_genomes_with_tag
taxnames_data = {} # taxname -> data {level (e.g. genus), taxid(e.g. adenovirus -> 1337), num_datasets}
if parse_metadata_from_ncbi and 'parse from entrez':
    print('Downloading taxonomy metadata from NCBI',flush=True)
    # Get taxdata from ncbi (in function)
    taxid_metadata_ncbi,taxids_missing_data = get_taxon_data_from_entrez(taxon_ids_all)
    #/
    # Determine name-variants at each taxnomic level
    taxids_handled = set() # IDE: check which taxids were ever assigned
    for taxid,data in sorted(taxid_metadata_ncbi.items()):
        if taxid in taxids_missing_data: continue # skip taxids with missing data
        for taxlevel,taxlevel_data in data['taxdata'].items():
            taxname = taxlevel_data['name']
            # save dataset as "present" to all its taxlevels
            if not taxlevel in taxid_metadata_ncbi_namesPerTaxLevel:            taxid_metadata_ncbi_namesPerTaxLevel[taxlevel] = {}
            if not taxname in taxid_metadata_ncbi_namesPerTaxLevel[taxlevel]:   taxid_metadata_ncbi_namesPerTaxLevel[taxlevel][taxname] = 0
            taxid_metadata_ncbi_namesPerTaxLevel[taxlevel][taxname] += 1
            #/
            # For defined taxlevelnames (i.e. when not "None") save data at taxname (e.g. adenovirus->{taxid,taxlevel,etc})
            if taxname != None:
                taxlevel_taxid = taxlevel_data['taxid']
                # bugcheck if this taxname as already added but does not have the same taxid as before. This is not expected...
                if taxname in taxnames_data and taxnames_data[taxname]['taxid'] != taxlevel_taxid:
                    print('FATAL: Already saved a taxid to this taxname and taxlevel. The new value does not match the previous!')
                    print('Previous value:'+ taxname+'='+str(taxnames_data[taxname]['taxid']))
                    print('New value:'+ taxname+'='+str(taxlevel_taxid))
                    print('Terminating!')
                    sys.exit()
                #/
                # save
                taxids_handled.add(taxid)
                if not taxname in taxnames_data:
                    taxnames_data[taxname] = {'taxid':taxlevel_taxid,'level':taxlevel,'datasets':set()}
                #/
            #/
        # IDE: check if input ever has "None" in species and serotype and "no rank"
        keys_to_check = ['species','serotype','no rank']
        keys_with_values_taxids = []
        for key_to_check in keys_to_check:
            if data['taxdata'][key_to_check]['taxid'] != None:
                keys_with_values_taxids.append(data['taxdata'][key_to_check]['taxid'])
        #/
        # Check if current taxid was ever found in any of parsed taxlevel columns. If not, then add custom "undefined_strain" taxlevelname
        taxid_matched = False
        for taxlevelname in data['taxdata']:
            if data['taxdata'][taxlevelname]['taxid'] != None and data['taxdata'][taxlevelname]['taxid'] == taxid:
                taxid_matched = True
        if not taxid_matched:
            print('FATAL: Taxid did not have a home!')
            print('Terminating!')
            sys.exit()
        #/
    #/
    #/
    # Add datasets in "taxnames_data"
    for taxid,data in sorted(taxid_metadata_ncbi.items()):
        for taxlevel,taxlevel_data in data['taxdata'].items():
            taxname = taxlevel_data['name']
            # skip if "none"
            if taxname == None: continue
            #/
            for bvbrc_id in taxon_ids_all_bvbrc_ids[taxid]:
                taxnames_data[taxname]['datasets'].add(bvbrc_id)
    #/
#/
##/
## Generate "BVBRC" metadata (info parsed from bvbrc summary/metadata files)
# check if user wants to download the bvbrc summary-file to parse additional data
if 0 and 'parse summary-file':
    asd = 1
#/
##/
##/

## Output metadata file "BASE" (info parsed from bvbrc fasta)
print('Generating output: metadata_base.tsv',flush=True)
with open(download_workdir+'/'+'metadata_base.tsv','w') as nf:
    # write header
    header = ['taxid','name','num','seq_min','seq_max','seq_med','seq_mean','description']
    nf.write('\t'.join(header)+'\n')
    #/
    # write data (sorted by number of datasets found which passed filtering criteria)
    for taxid,data in sorted(taxid_metadata_base.items(),key=lambda x: x[1]['num'], reverse=True):
        writeArr = [ taxid,data['desc']['species'],data['num'],data['min'],data['max'],data['median'],data['mean'],data['desc']['full'] ]
        nf.write('\t'.join(map(str,writeArr))+'\n')
    #/
##/
## Output metadata file "NCBI" (info parsed from NCBI via Entrez) PLUS number of genomes with+without tag-passed seqs
if parse_metadata_from_ncbi:
    print('Generating output: metadata_ncbi.tsv',flush=True)
    with open(download_workdir+'/'+'metadata_ncbi.tsv','w') as nf:
        # write header
        header = ['taxid','name','num_with_tag'] + taxnames_to_parse_ordered + [x+'_taxid' for x in taxnames_to_parse_ordered]
        nf.write('\t'.join(header)+'\n')
        #/
        # write data (sorted by number of genomes with "tag" [e.g. "complete genome"])
        for taxid,data in sorted(taxid_metadata_ncbi.items(),key=lambda x: 0 if not x[0] in taxid_metadata_base else taxid_metadata_base[x[0]]['num'],reverse=True):
            # get dataset counts for taxid from "BASE" metadata
            taxid_num_datasets = None
            if taxid in taxid_metadata_base:
                taxid_num_datasets = taxid_metadata_base[taxid]['num']
            #/
            ## compile writeArr
            # first columns
            writeArr = [taxid,data['name'],taxid_num_datasets]
            #/
            # taxon name columns
            for taxname in taxnames_to_parse_ordered:
                taxname_value = data['taxdata'][taxname]['name']
                writeArr.append(taxname_value)
            #/
            # taxon ID columns
            for taxname in taxnames_to_parse_ordered:
                taxname_value = data['taxdata'][taxname]['taxid']
                writeArr.append(taxname_value)
            nf.write('\t'.join(map(str,writeArr))+'\n')
            #/
            ##/
        #/
##/
## Output metadata-summary table: Columns per each taxlevel with 2 sub-columns: (1) taxlevelname (2) number of datasets with taxname
if parse_metadata_from_ncbi and 'parse from entrez':
    print('Generating output: metadata_ncbi_summary.tsv',flush=True)
    with open(download_workdir+'/'+'metadata_ncbi_summary.tsv','w') as nf:
        # write header
        header = []
        for taxlevelname in taxnames_to_parse_ordered:
            header.append(taxlevelname)
            header.append('taxid')
            header.append('count')
            header.append('count_with_tag')
        nf.write('\t'.join(header)+'\n')
        #/
        # compile output columns (will transpose to rows later)
        cols = []
        for taxlevelname in taxnames_to_parse_ordered:
            cols.append([]) # init name-col
            cols.append([]) # init taxid-col
            cols.append([]) # init counts-col
            cols.append([]) # init counts_with_tag-col
            
            if not taxlevelname in taxid_metadata_ncbi_namesPerTaxLevel: continue # skip if no data at this taxonomic level
            
            # add column rows (sorted by counts)
            for taxname,count in sorted(taxid_metadata_ncbi_namesPerTaxLevel[taxlevelname].items(),key=lambda x: x[1],reverse=True):
                # if none then skip
                if taxname == None: continue
                #/
                bvbrc_ids_num = len(taxnames_data[taxname]['datasets'])
                cols[-4].append(taxname)
                cols[-3].append(taxnames_data[taxname]['taxid'])
                cols[-2].append(bvbrc_ids_num)
                cols[-1].append(count)
            
            if 0 and 'old number-based saver':
                for taxname,count in sorted(taxid_metadata_ncbi_namesPerTaxLevel[taxlevelname].items(),key=lambda x: x[1],reverse=True):
                    cols[-2].append(taxname)
                    cols[-1].append(count)
            #/
        #/
        # write columns as rows
        rowidx = 0
        while True:
            # get row values from columns
            row = []
            for col in cols:
                if len(col) > rowidx:
                    row.append(col[rowidx])
                else:
                    row.append('')
            #/
            # check if row is all empty, that means we are finished
            if row.count('') == len(row):
                break
            #/
            # else, write to file and iterate row
            nf.write('\t'.join(map(str,row))+'\n')
            rowidx += 1
            #/
        #/
##/
## Dump taxname -> bvbrc dict
if parse_metadata_from_ncbi and 'parse from entrez':
    print('Generating output: metadata_ncbi_taxnames_data.pickle',flush=True)
    with open(download_workdir+'/'+'metadata_ncbi_taxnames_data.pickle','wb') as nf:
        pickle.dump(taxnames_data,nf)
##/
###/

### Copy-out metadata files for user to outer workdir
print('Copy-out metadata files to workdir')
files_to_copy = ['metadata_base.tsv','metadata_ncbi.tsv','metadata_ncbi_summary.tsv']
for file_to_copy in files_to_copy:
    # copy-out files with prepended "virus_"
    shutil.copy2(download_workdir+'/'+file_to_copy,workdir+'/'+'virus_'+file_to_copy)
###/