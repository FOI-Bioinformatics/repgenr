#!/usr/bin/env python3

import os
import sys
import urllib.request
import tarfile
import argparse
from matplotlib import pyplot as plt

taxonomy_ordered = ('domain','phylum','class','family','genus','species')

### Parse input arguments
# setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-r','--release',required=True,help='Release-version of GTDB to use, must include major followed by minor release versions.\nExample: 207.0, where major release is 207 and minor release is 0\nSee available releases at https://data.gtdb.ecogenomic.org/releases/')
parser.add_argument('-v','--version',required=True,choices=('bac120','ar53',),help='Database version to use')
parser.add_argument('-d','--dataset',required=True,choices=('all','rep',),help='Dataset to use')

parser.add_argument('-l','--level',required=True,choices=('family','genus','species',),help='Taxonomic level to fetch organisms from (e.g. family) related to target')
parser.add_argument('-tf','--target_family',help='Target family (Example: francisellaceae)')
parser.add_argument('-tg','--target_genus',help='Target genus (Example: francisella)')
parser.add_argument('-ts','--target_species',help='Target species (Example: tularensis)')

parser.add_argument('-wd','--workdir',required=True,help='Path to working directory. This folder will be created if not already present')

parser.add_argument('--nodownload',help='Use this flag if the database has already been downloaded (to prevent unneccessary load to GTDB)',action='store_true')

parser.add_argument('--outgroup_accession',default=None,help='Provide an NCBI accession number to use as outgroup in taxonomic output (e.g. GCA_01234567.8')

#/
# parse input
args = parser.parse_args()
#args = parser.parse_args(['-r', '207.0', '-v', 'bac120', '-d', 'all', '-l', 'family', '-tg', 'francisella', '-ts', 'tularensis', '-wd', 'tularensis2', '--nodownload'])

gtdb_release = args.release
gtdb_version = args.version
gtdb_dataset = args.dataset

level = args.level.lower()
target_family = args.target_family
target_genus = args.target_genus
target_species = args.target_species

workdir = args.workdir

nodownload = args.nodownload

outgroup_accession = args.outgroup_accession
#/
# validate input
if gtdb_release.find('.') == -1:
    print('Could not parse major and minor release! Expected format: 207.0, where major release is 207 and minor release is 0.')
    print('Input value: '+str(gtdb_release))
    sys.exit()
    
try:
    gtdb_release = float(gtdb_release)
except:
    print('Could not parse major and minor release! Expected format: 207.0, where major release is 207 and minor release is 0.')
    print('Input value: '+str(gtdb_release))
    sys.exit()

if not (target_genus or target_family):
    print('Must supply a genus or family as input using --target_genus or --target_family (e.g., francisella or francisellaceae)')
    sys.exit()
    
if target_genus and target_family:
    if not target_species:
        print('Both family and genus was supplied as input. Please supply either genus or family as input using --target_genus or --target_family (e.g., francisella or francisellaceae)')
    else:
        print('Both family and genus was supplied as input. Please supply only genus in addition to --target_species as input using --target_genus (e.g., francisella)')
    sys.exit()

if target_species and not target_genus:
    print('Must supply a target genus for species using --target_genus (e.g., francisella)')
    sys.exit()
#/
# Final formatting
if target_family:           target_family = target_family.lower()
if target_genus:            target_genus = target_genus.lower()
if target_species:          target_species = target_species.lower()
#/
###/

### Download gtdb database
print('Downloading metadata-file from GTDB...')
if not os.path.exists(workdir):     os.makedirs(workdir)

metadata_target = 'https://data.gtdb.ecogenomic.org/releases/release'+str(int(gtdb_release))+'/'+str(gtdb_release)+'/'+gtdb_version+'_metadata_r'+str(int(gtdb_release))+'.tar.gz'
metadata_file = os.path.basename(metadata_target)
try:
    if not nodownload:
        urllib.request.urlretrieve(metadata_target, workdir+'/'+metadata_file)
    else:
        print(' ^ Will use downloaded database, flag --nodownload specified')
except:
    print('Could not fetch metadata from:')
    print(metadata_target)
    print('Make sure that the release version, database version, and dataset version exists.')
    sys.exit()
###/

### Read GTDB database
print('Parsing metadata-file...')

fo_tar = tarfile.open(workdir+'/'+metadata_file)
accessions_data = {}
for content in fo_tar.getmembers():
    content_handle = fo_tar.extractfile(content)

    # We expect a .tsv file    
    #if content_handle.name.find('.tsv') == -1: continue
    #/
    
    header = None
    for enum_file,line_raw in enumerate(content_handle):
        line = line_raw.decode('utf-8')
        line = line.strip('\n')
        line = line.split('\t')
        # parse header
        if enum_file == 0:
            header = line
            continue
        #/
        ## parse lines
        accession = None
        accession_ncbi = None
        genome_representative_gtdb = None
        representative_gtdb = None
        
        tax_gtdb = {}
        tax_ncbi = {}
        for tax in taxonomy_ordered:
            tax_gtdb[tax] = None
            tax_ncbi[tax] = None
            
        strain_ID_ncbi = ''
        for enum,entry in enumerate(line):
            column = header[enum]
            
            # parse accession
            if column == 'accession':
                accession = entry.replace('GB_','').replace('RS_','')
            if column == 'ncbi_genbank_assembly_accession':
                accession_ncbi = entry
            #/
            # parse ...?
            if column == 'gtdb_genome_representative':
                genome_representative_gtdb = entry
            if column == 'gtdb_representative':
                representative_gtdb = entry
            #/
            # parse gtdb taxonomy
            if column == 'gtdb_taxonomy':
                entry = entry.split(';')
                for chunk in entry:
                    for tax in tax_gtdb:
                        check_key = tax[0]+'__'
                        
                        if chunk.find(check_key) != -1:
                            tax_gtdb[tax] = chunk.replace(check_key,'')
            #/
            # parse ncbi taxonomy
            if column == 'ncbi_taxonomy':
                entry = entry.split(';')
                for chunk in entry:
                    for tax in tax_ncbi:
                        check_key = tax[0]+'__'
                        
                        if chunk.find(check_key) != -1:
                            tax_ncbi[tax] = chunk.replace(check_key,'')
            #/
            # parse ncbi strain identifiers
            if column == 'ncbi_strain_identifiers':
                if entry != 'none':
                    strain_ID_ncbi = entry
            #/
  
        # bugcheck if we have multiple entries for accession
        if accession in accessions_data:
            print('Found another entry for accession?! I expect only one ',accession)
            sys.exit()
        #/
        # compile save data
        save_data = {'accession':accession,'accession_ncbi':accession_ncbi,
                     'tax_gtdb':tax_gtdb,'tax_ncbi':tax_ncbi,
                     'strain_ID_ncbi':strain_ID_ncbi}
        #/
        # remove spaces and redundant information in "species" taxonomy tag
        for db in ('gtdb','ncbi',):
            save_data['tax_'+db]['species'] = save_data['tax_'+db]['species'].replace(save_data['tax_'+db]['genus'],'')
            save_data['tax_'+db]['species'] = save_data['tax_'+db]['species'].replace(' ','')
            save_data['tax_'+db]['species'] = save_data['tax_'+db]['species'].replace('_','-')
        #/
        # save to outer
        accessions_data[accession] = save_data
        #/
        ##/
fo_tar.close()

# check if accessions are always GCA/GCF
accessions_prefixes = set()
for acc in accessions_data:
    accessions_prefixes.add(acc.split('_')[0])
#/
###/

### Parse specified organism and get all other samples at specified level
print('Parsing target and other organisms at specified level...')
## Get target values at taxonomic levels
target_value_at_level = None
for acc,data in accessions_data.items():
    # Check if accession matches our target
    #@ check if we have species + genus input
    target_found = False
    if target_species and target_genus:
        if target_species.lower() == data['tax_gtdb']['species'].lower() and target_genus.lower() == data['tax_gtdb']['genus'].lower():
            target_found = True
    #@ else do without species
    if (target_genus or target_family) and not target_species:
        if target_genus and target_genus.lower() == data['tax_gtdb']['genus'].lower():
            target_found = True
        if target_family and target_family.lower() == data['tax_gtdb']['family'].lower():
            target_found = True
    #/
    # if match to target, get values at all specified
    if target_found:
        target_value_at_level = data['tax_gtdb'][level]
    #/
    
if not target_value_at_level:
    print('Did not find any entries in database for input target family="'+str(target_family)+'" genus="'+str(target_genus)+'" species="'+str(target_species)+'" at taxonomic level "'+level+'"')
    print('Please ensure that your target exists in the database!')
    sys.exit()
##/
## Get other organisms at specified level
accessions_selected = {}
for acc,data in accessions_data.items():
    # check if taxonomic level matches the target value at the same level
    levels_match = False
    if data['tax_gtdb'][level] == target_value_at_level:
        accessions_selected[acc] = data
    #/
##/
## Get organism at upper taxonomic level (for rooting phylogenetic tree by using an outgroup organism)
tax_level_upper_accessions = {}
# Parse manually assigned accession
if outgroup_accession:
    tax_level_upper_accessions[outgroup_accession] = True
    
else: # else, get one via metadata
    # Find upper level
    tax_level_upper = None
    for i in range(len(taxonomy_ordered)):
        if i < len(taxonomy_ordered):
            lower_tax = taxonomy_ordered[i+1]
            upper_tax = taxonomy_ordered[i]
            
            if lower_tax == level:
                tax_level_upper = upper_tax
                break
    #/
    # Determine what the upper level value is in the target
    tax_level_upper_val = None
    for any_acc,any_data in accessions_selected.items():
        tax_level_upper_val = any_data['tax_gtdb'][tax_level_upper]
        break
    #/
    # Get a sample at val_level_upper that is not part of target set
    for any_acc,any_data in accessions_data.items():
        if any_acc in accessions_selected: continue # skip if part of selected targets
        if any_data['tax_gtdb'][tax_level_upper] == tax_level_upper_val:
            tax_level_upper_accessions[any_acc] = any_data
            
            # Only get X amount of higher-tax organisms
            if len(tax_level_upper_accessions) >= 1:
                break
            #/
    #/
##/
###/

### Summarize selected metadata
print('Summarizing data...')
tax_val_counts = {}
for acc, data in accessions_selected.items():
    for tax,val in data['tax_gtdb'].items():
        if not tax in tax_val_counts:       tax_val_counts[tax] = {}
        if not val in tax_val_counts[tax]:  tax_val_counts[tax][val] = 0
        tax_val_counts[tax][val] += 1
###/

### Output summary and selected accessions
print('Writing output...')
if not os.path.exists(workdir):          os.makedirs(workdir)
## Summary
# table
tax_order = ('family','genus','species',)
yvals = []
xvals = []
with open(workdir+'/'+'metadata_summary.tsv','w') as nf:
    header = ['level','name','abundance']
    nf.write('\t'.join(map(str,header))+'\n')
    for tax in tax_order:
        yvals.append([])
        xvals.append([])
        for val,abu in sorted(tax_val_counts[tax].items(),key=lambda x: x[1], reverse=True):
            yvals[-1].append(abu)
            xvals[-1].append(val)
            writeArr = [tax,val,abu]
            nf.write('\t'.join(map(str,writeArr))+'\n')
#/
# barplot, number in each level
x_vals = xvals[0] + xvals[1] + xvals[2]
fig, ax = plt.subplots()
bar_fam = ax.bar(xvals[0],yvals[0],label='family')
bar_genus = ax.bar(xvals[1],yvals[1],label='genus')
bar_species = ax.bar(xvals[2],yvals[2],label='species')
ax.xaxis.set_ticks(list(range(len(x_vals))))
ax.xaxis.set_ticklabels(x_vals,rotation=90)
plt.legend(loc='upper right')
plt.title('Number of datasets in each taxonomic level')
plt.tight_layout()
plt.savefig(workdir+'/'+'metadata_summary_number_in_level.png',dpi=200)
#/
# barplot, number at each level
yvals2 = list(map(len,yvals))
xvals2 = ['family','genus','species']
fig, ax = plt.subplots()
bar = ax.bar(xvals2,yvals2)
ax.bar_label(bar)
ax.xaxis.set_ticks(list(range(len(xvals2))))
ax.xaxis.set_ticklabels(xvals2,rotation=90)
plt.title('Number of taxa per taxonomic level')
plt.tight_layout()
plt.savefig(workdir+'/'+'metadata_summary_number_per_level.png',dpi=200)
#/
##/

## Metadata for selected accessions
with open(workdir+'/'+'metadata_selected.tsv','w') as nf:
    nf.write(str(accessions_selected))
##/

## Accession numbers for selected accessions
if 0 :
    with open(workdir+'/'+'accessions_selected.tsv','w') as nf:
        nf.write('\n'.join(sorted(accessions_selected)))
##/

## Output level
with open(workdir+'/'+'metadata_level.txt','w') as nf:
    nf.write(level+'\n')
##/

## Output outgroup organism
with open(workdir+'/'+'outgroup_accession.txt','w') as nf:
    nf.write(list(tax_level_upper_accessions)[0]+'\n')
##/
###/