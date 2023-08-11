#!/usr/bin/env python3

import os
import sys
import argparse


def summarize_derep_genomes(workdir):
    # Get genome representants
    genome_representants = []
    for file_ in os.listdir(workdir+'/'+'genomes_derep_representants'):
        if file_.endswith('.fasta') or file_.endswith('.fasta.gz'):
            genome_representants.append(file_)
    #/
    # Get dereplicated-and-removed genomes from "stage 2"/chunks (if applicable)
    dereps_stage2 = {}
    if os.path.exists(workdir+'/'+'derep_chunks_clustered_genomes.tsv'):
        with open(workdir+'/'+'derep_chunks_clustered_genomes.tsv','r') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split('\t')
                derep_parent,derep_child,chunk = line
                if not derep_parent in dereps_stage2:       dereps_stage2[derep_parent] = set()
                dereps_stage2[derep_parent].add(derep_child)
    #/
    # Get dereplicated-and-removed genomes from "stage 1"
    dereps_stage1 = {}
    with open(workdir+'/'+'derep_clustered_genomes.tsv','r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split('\t')
            derep_parent,derep_child = line
            if not derep_parent in dereps_stage1:       dereps_stage1[derep_parent] = set()
            dereps_stage1[derep_parent].add(derep_child)
    #/
    # Summarize all children under every genome representant
    dereps_summarized = {}
    for rep in genome_representants:
        dereps_summarized[rep] = set()
        
        # check if this rep has childs in stage1
        if rep in dereps_stage1:
            stage1_childs = dereps_stage1[rep]
            dereps_summarized[rep].update(stage1_childs)
            # check if any stage1-childs had stage2-childs
            for stage1_child in stage1_childs:
                if stage1_child in dereps_stage2:
                    stage2_childs = dereps_stage2[stage1_child]
                    dereps_summarized[rep].update(stage2_childs)
            #/
        #/
        # Check if this rep had childs in stage2
        if rep in dereps_stage2:
            stage2_childs = dereps_stage2[rep]
            dereps_summarized[rep].update(stage2_childs)
        #/
    #/
    # BUGCHECK: Confirm that every child only has one parent
    derep_childs_parents = {}
    for parent,childs in dereps_summarized.items():
        for child in childs:
            if not child in derep_childs_parents:       derep_childs_parents[child] = set()
            derep_childs_parents[child].add(parent)
        
    for child,parent in derep_childs_parents.items():
        if len(parent) > 1:
            sys.exit('FATAL: There was a bug: a child had multiple parents during summary of dereplication!')
    #/
    # Write output
    with open(workdir+'/'+'derep_genomes_summary.tsv','w') as nf:
        # write header
        header = ['representant','num_contained','contained']
        nf.write('\t'.join(map(str,header))+'\n')
        #/
        # write rows
        for rep,childs in dereps_summarized.items():
            writeArr = [rep,len(childs),'\t'.join(sorted(childs))]
            nf.write('\t'.join(map(str,writeArr))+'\n')
        #/
    #/
    # Finalize
    print('Dereplication summary written to "derep_genomes_summary.tsv"')
    #/

####/

### Execute code o this script if it was called stand-alone
if __name__ == "__main__":
    # setup
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-wd','--workdir',required=True,help='Path to working directory')
    #/
    # parse input
    args = parser.parse_args()
    
    workdir = args.workdir
    #/
    # validate input
    if not os.path.exists(workdir):
        print('Could not locate working directory at:')
        print(workdir)
        sys.exit()
    
    summarize_derep_genomes(workdir)
###/