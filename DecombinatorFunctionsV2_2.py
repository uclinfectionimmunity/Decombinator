#!/usr/bin/python

print 'Loading...'

import sys, argparse, os
import numpy as np
import decimal as dec
import string
import operator as op
import collections as coll
from Bio import SeqIO
import time
from string import Template
from operator import itemgetter, attrgetter
import Levenshtein as lev
from acora import AcoraBuilder
import platform

v_half_split, j_half_split = [10,6] # Do not change - V tags are split at position 10, J at position 10, to look for half tags if no full tag is found.

def create_folder(outputfile):

    currentpath = os.getcwd()
    
    if platform.system() == 'Windows':
        newpath = currentpath+'\\results_'+str(outputfile)+'\\' ## Ensure correct for specified platform
        if not os.path.exists(newpath):
            os.makedirs(newpath)
    elif platform.system() == 'Linux':
        newpath = currentpath+'/results_'+str(outputfile)+'/' ## Ensure correct for specified platform
        if not os.path.exists(newpath):
            os.makedirs(newpath)
    elif platform.system() == 'Darwin':
        newpath = currentpath+'/results_'+str(outputfile)+'/' ## Ensure correct for specified platform
        if not os.path.exists(newpath):
            os.makedirs(newpath)
            
    return newpath

def analysis( inputfile, outputfile, with_reverse_complement_search, species, chain, barcode, barcodestart1, barcodeend1, barcodestart2, barcodeend2, newpath, omitN=True ):
    
    out_file_alpha = open( newpath+outputfile+'_alpha.txt',"w")
    out_file_beta = open( newpath+outputfile+'_beta.txt',"w")
    out_file_gamma = open( newpath+outputfile+'_gamma.txt',"w")
    out_file_delta = open( newpath+outputfile+'_delta.txt',"w")
    log_file = open( newpath+outputfile+'_summary.txt',"w")
    Nseqs = 0

    # Import all known V and J alpha, beta, gamma and delta gene sequences from IMGT
    ################
    
    if species == 'human':
        va_regions, vb_regions, vg_regions, vd_regions, ja_regions, jb_regions, jg_regions, jd_regions = load_human_gene_sequences()
        
        va_seqs, half1_va_seqs, half2_va_seqs, jump_to_end_va = get_v_tags(open("humantags_trav.txt", "rU"), v_half_split)
        vb_seqs, half1_vb_seqs, half2_vb_seqs, jump_to_end_vb = get_v_tags(open("humantags_trbv.txt", "rU"), v_half_split)
        vg_seqs, half1_vg_seqs, half2_vg_seqs, jump_to_end_vg = get_v_tags(open("humantags_trgv.txt", "rU"), v_half_split)
        vd_seqs, half1_vd_seqs, half2_vd_seqs, jump_to_end_vd = get_v_tags(open("humantags_trdv.txt", "rU"), v_half_split)
        
        ja_seqs, half1_ja_seqs, half2_ja_seqs, jump_to_start_ja = get_j_tags(open("humantags_traj.txt", "rU"), j_half_split)
        jb_seqs, half1_jb_seqs, half2_jb_seqs, jump_to_start_jb = get_j_tags(open("humantags_trbj.txt", "rU"), j_half_split)
        jg_seqs, half1_jg_seqs, half2_jg_seqs, jump_to_start_jg = get_j_tags(open("humantags_trgj.txt", "rU"), j_half_split)
        jd_seqs, half1_jd_seqs, half2_jd_seqs, jump_to_start_jd = get_j_tags(open("humantags_trdj.txt", "rU"), j_half_split)
        
    elif species == 'mouse':
        va_regions, vb_regions, vg_regions, vd_regions, ja_regions, jb_regions, jg_regions, jd_regions = load_mouse_gene_sequences()
        
        va_seqs, half1_va_seqs, half2_va_seqs, jump_to_end_va = get_v_tags(open("mousetags_trav.txt", "rU"), v_half_split)
        vb_seqs, half1_vb_seqs, half2_vb_seqs, jump_to_end_vb = get_v_tags(open("mousetags_trbv.txt", "rU"), v_half_split)
        vg_seqs, half1_vg_seqs, half2_vg_seqs, jump_to_end_vg = get_v_tags(open("mousetags_trgv.txt", "rU"), v_half_split)
        vd_seqs, half1_vd_seqs, half2_vd_seqs, jump_to_end_vd = get_v_tags(open("mousetags_trdv.txt", "rU"), v_half_split)
        
        ja_seqs, half1_ja_seqs, half2_ja_seqs, jump_to_start_ja = get_j_tags(open("mousetags_traj.txt", "rU"), j_half_split)
        jb_seqs, half1_jb_seqs, half2_jb_seqs, jump_to_start_jb = get_j_tags(open("mousetags_trbj.txt", "rU"), j_half_split)
        jg_seqs, half1_jg_seqs, half2_jg_seqs, jump_to_start_jg = get_j_tags(open("mousetags_trgj.txt", "rU"), j_half_split)
        jd_seqs, half1_jd_seqs, half2_jd_seqs, jump_to_start_jd = get_j_tags(open("mousetags_trdj.txt", "rU"), j_half_split)

    ### Build keyword tries using V and J tags for fast assignment
    
    va_key = build_keyword_tries(va_seqs)
    vb_key = build_keyword_tries(vb_seqs)
    vg_key = build_keyword_tries(vg_seqs)
    vd_key = build_keyword_tries(vd_seqs)

    ja_key = build_keyword_tries(ja_seqs)
    jb_key = build_keyword_tries(jb_seqs)
    jg_key = build_keyword_tries(jg_seqs)
    jd_key = build_keyword_tries(jd_seqs)    

    ### Build keyword tries for first and second halves of both V and J tags
    
    half1_va_key = build_keyword_tries(half1_va_seqs)
    half1_vb_key = build_keyword_tries(half1_vb_seqs)
    half1_vg_key = build_keyword_tries(half1_vg_seqs)
    half1_vd_key = build_keyword_tries(half1_vd_seqs)

    half1_ja_key = build_keyword_tries(half1_ja_seqs)
    half1_jb_key = build_keyword_tries(half1_jb_seqs)
    half1_jg_key = build_keyword_tries(half1_jg_seqs)
    half1_jd_key = build_keyword_tries(half1_jd_seqs)

    half2_va_key = build_keyword_tries(half2_va_seqs)
    half2_vb_key = build_keyword_tries(half2_vb_seqs)
    half2_vg_key = build_keyword_tries(half2_vg_seqs)
    half2_vd_key = build_keyword_tries(half2_vd_seqs)

    half2_ja_key = build_keyword_tries(half2_ja_seqs)
    half2_jb_key = build_keyword_tries(half2_jb_seqs)
    half2_jg_key = build_keyword_tries(half2_jg_seqs)
    half2_jd_key = build_keyword_tries(half2_jd_seqs)

    ### Initialise variables
    
    assigned_count_alpha = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
    assigned_count_beta = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
    assigned_count_gamma = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
    assigned_count_delta = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
    seq_count = 0 # this will simply track the number of sequences analysed in file
    error0_count = 0 # number of sequences with no errors in V tags
    error1_count = 0 # number of sequences with 1 error in V tags
    t0 = time.time() # Begin timer

    ### Begin analysing sequences
    
    item = str(inputfile)
    print 'Importing sequences from', item,'and assigning V and J regions...'
    handle = open(item, "rU")
              
    for record in SeqIO.parse(handle, "fastq"):
        seq_count += 1
        print seq_count
        found_seq_match = 0
                  
        ## DETERMINE BARCODE SEQUENCE AT START OF SEQUENCE
        if barcode == True:
            barcode_seq = str(record.seq)[barcodestart1:barcodeend1]+str(record.seq)[barcodestart2:barcodeend2]
            barcode_qual = str(record.format("fastq")).split("\n")[3][barcodestart1:barcodeend1]+str(record.format("fastq")).split("\n")[3][barcodestart2:barcodeend2]
            
        else:
            barcode_seq = 'NA'
            barcode_qual = 'NA'
            
        if 'N' in str(record.seq):
            Nseqs += 1
            
        if ((chain == "alpha") or (chain == "all")) and found_seq_match==0:
            assigned_count_alpha, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record.seq), str(record.id),
                                                                                                        seq_count, assigned_count_alpha, Nseqs,
                                                                                                        va_key, ja_key,
                                                                                                        va_seqs, ja_seqs,
                                                                                                        half1_va_seqs, half2_va_seqs,
                                                                                                        half1_ja_seqs, half2_ja_seqs,
                                                                                                        jump_to_end_va, jump_to_start_ja,
                                                                                                        va_regions, ja_regions,
                                                                                                        half1_va_key, half2_va_key,
                                                                                                        half2_ja_key, half2_ja_key,
                                                                                                        out_file_alpha,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
        if ((chain == "beta") or (chain == "all")) and found_seq_match==0:
            assigned_count_beta, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record.seq), str(record.id),
                                                                                                        seq_count, assigned_count_beta, Nseqs,
                                                                                                        vb_key, jb_key,
                                                                                                        vb_seqs, jb_seqs,
                                                                                                        half1_vb_seqs, half2_vb_seqs,
                                                                                                        half1_jb_seqs, half2_jb_seqs,
                                                                                                        jump_to_end_vb, jump_to_start_jb,
                                                                                                        vb_regions, jb_regions,
                                                                                                        half1_vb_key, half2_vb_key,
                                                                                                        half2_jb_key, half2_jb_key,
                                                                                                        out_file_beta,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
        if ((chain == "gamma") or (chain == "all")) and found_seq_match==0:
            assigned_count_gamma, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record.seq), str(record.id),
                                                                                                        seq_count, assigned_count_gamma, Nseqs,
                                                                                                        vg_key, jg_key,
                                                                                                        vg_seqs, jg_seqs,
                                                                                                        half1_vg_seqs, half2_vg_seqs,
                                                                                                        half1_jg_seqs, half2_jg_seqs,
                                                                                                        jump_to_end_vg, jump_to_start_jg,
                                                                                                        vg_regions, jg_regions,
                                                                                                        half1_vg_key, half2_vg_key,
                                                                                                        half2_jg_key, half2_jg_key,
                                                                                                        out_file_gamma,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
        if ((chain == "delta") or (chain == "all")) and found_seq_match==0:
            assigned_count_delta, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record.seq), str(record.id),
                                                                                                        seq_count, assigned_count_delta, Nseqs,
                                                                                                        vd_key, jd_key,
                                                                                                        vd_seqs, jd_seqs,
                                                                                                        half1_vd_seqs, half2_vd_seqs,
                                                                                                        half1_jd_seqs, half2_jd_seqs,
                                                                                                        jump_to_end_vd, jump_to_start_jd,
                                                                                                        vd_regions, jd_regions,
                                                                                                        half1_vd_key, half2_vd_key,
                                                                                                        half2_jd_key, half2_jd_key,
                                                                                                        out_file_delta,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
        if found_seq_match == 0 and with_reverse_complement_search == True:
            
            #####################
            ## REVERSE COMPLEMENT
            #####################
            
            record_reverse = record.reverse_complement()
            
            if ((chain == "alpha") or (chain == "all")) and found_seq_match==0:
                assigned_count_alpha, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record_reverse.seq), str(record.id),
                                                                                                        seq_count, assigned_count_alpha, Nseqs,
                                                                                                        va_key, ja_key,
                                                                                                        va_seqs, ja_seqs,
                                                                                                        half1_va_seqs, half2_va_seqs,
                                                                                                        half1_ja_seqs, half2_ja_seqs,
                                                                                                        jump_to_end_va, jump_to_start_ja,
                                                                                                        va_regions, ja_regions,
                                                                                                        half1_va_key, half2_va_key,
                                                                                                        half2_ja_key, half2_ja_key,
                                                                                                        out_file_alpha,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
                                
            if ((chain == "beta") or (chain == "all")) and found_seq_match==0:
                assigned_count_beta, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record_reverse.seq), str(record.id),
                                                                                                        seq_count, assigned_count_beta, Nseqs,
                                                                                                        vb_key, jb_key,
                                                                                                        vb_seqs, jb_seqs,
                                                                                                        half1_vb_seqs, half2_vb_seqs,
                                                                                                        half1_jb_seqs, half2_jb_seqs,
                                                                                                        jump_to_end_vb, jump_to_start_jb,
                                                                                                        vb_regions, jb_regions,
                                                                                                        half1_vb_key, half2_vb_key,
                                                                                                        half2_jb_key, half2_jb_key,
                                                                                                        out_file_beta,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
                                
            if ((chain == "gamma") or (chain == "all")) and found_seq_match==0:
                assigned_count_gamma, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record_reverse.seq), str(record.id),
                                                                                                        seq_count, assigned_count_gamma, Nseqs,
                                                                                                        vg_key, jg_key,
                                                                                                        vg_seqs, jg_seqs,
                                                                                                        half1_vg_seqs, half2_vg_seqs,
                                                                                                        half1_jg_seqs, half2_jg_seqs,
                                                                                                        jump_to_end_vg, jump_to_start_jg,
                                                                                                        vg_regions, jg_regions,
                                                                                                        half1_vg_key, half2_vg_key,
                                                                                                        half2_jg_key, half2_jg_key,
                                                                                                        out_file_gamma,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
            
            if ((chain == "delta") or (chain == "all")) and found_seq_match==0:
                assigned_count_delta, seq_count, Nseqs, found_seq_match, error0_count, error1_count  = engine(str(record_reverse.seq), str(record.id),
                                                                                                        seq_count, assigned_count_delta, Nseqs,
                                                                                                        vd_key, jd_key,
                                                                                                        vd_seqs, jd_seqs,
                                                                                                        half1_vd_seqs, half2_vd_seqs,
                                                                                                        half1_jd_seqs, half2_jd_seqs,
                                                                                                        jump_to_end_vd, jump_to_start_jd,
                                                                                                        vd_regions, jd_regions,
                                                                                                        half1_vd_key, half2_vd_key,
                                                                                                        half2_jd_key, half2_jd_key,
                                                                                                        out_file_delta,
                                                                                                        error0_count, error1_count,
                                                                                                        barcode_seq, barcode_qual
                                                                                                        )
            
    handle.close()
    out_file_alpha.close()
    out_file_beta.close()
    out_file_gamma.close()
    out_file_delta.close()

    ### Print analysis summary to log file

    timed = time.time() - t0
    print 'Completed analysis of sequences'
    print >> log_file, seq_count, 'sequences were analysed'
    print >> log_file, assigned_count_alpha, 'TcR alpha sequences were successfully assigned'
    print >> log_file, assigned_count_beta, 'TcR beta sequences were successfully assigned'
    print >> log_file, assigned_count_gamma, 'TcR gamma sequences were successfully assigned'
    print >> log_file, assigned_count_delta, 'TcR delta sequences were successfully assigned'
    if error0_count+error1_count == 0:
        print "No TcR sequences found!"
    else:
        print >> log_file, 1-(20*error0_count+19*error1_count)/float(20*(error0_count+error1_count)), 'upper bound on sequencing error rate'
    print >> log_file, Nseqs, 'sequences contained ambiguous N nucleotides'
    print >> log_file, 'Time taken =', timed, 'seconds'

### Functions

def load_human_gene_sequences():
    
    print ('Importing known V and J gene segments and tags...')

    handle = open("human_TRAV_region.fasta", "rU")
    va_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRBV_region.fasta", "rU")
    vb_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRGV_region.fasta", "rU")
    vg_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRDV_region.fasta", "rU")
    vd_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRAJ_region.fasta", "rU")
    ja_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRBJ_region.fasta", "rU")
    jb_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRGJ_region.fasta", "rU")
    jg_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRDJ_region.fasta", "rU")
    jd_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    va_regions = []
    for j in range(0, len(va_genes)):
        va_regions.append(string.upper(va_genes[j].seq))

    vb_regions = []
    for j in range(0, len(vb_genes)):
        vb_regions.append(string.upper(vb_genes[j].seq))

    vg_regions = []
    for j in range(0, len(vg_genes)):
        vg_regions.append(string.upper(vg_genes[j].seq))

    vd_regions = []
    for j in range(0, len(vd_genes)):
        vd_regions.append(string.upper(vd_genes[j].seq))
        
    ja_regions = []
    for j in range(0, len(ja_genes)):
        ja_regions.append(string.upper(ja_genes[j].seq))

    jb_regions = []
    for j in range(0, len(jb_genes)):
        jb_regions.append(string.upper(jb_genes[j].seq))

    jg_regions = []
    for j in range(0, len(jg_genes)):
        jg_regions.append(string.upper(jg_genes[j].seq))

    jd_regions = []
    for j in range(0, len(jd_genes)):
        jd_regions.append(string.upper(jd_genes[j].seq))

    return va_regions, vb_regions, vg_regions, vd_regions, ja_regions, jb_regions, jg_regions, jd_regions

def load_mouse_gene_sequences():
    
    print ('Importing known V and J gene segments and tags...')

    handle = open("mouse_TRAV_region.fasta", "rU")
    va_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRBV_region.fasta", "rU")
    vb_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRGV_region.fasta", "rU")
    vg_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRDV_region.fasta", "rU")
    vd_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRAJ_region.fasta", "rU")
    ja_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRBJ_region.fasta", "rU")
    jb_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRGJ_region.fasta", "rU")
    jg_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("mouse_TRDJ_region.fasta", "rU")
    jd_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    va_regions = []
    for j in range(0, len(va_genes)):
        va_regions.append(string.upper(va_genes[j].seq))

    vb_regions = []
    for j in range(0, len(vb_genes)):
        vb_regions.append(string.upper(vb_genes[j].seq))

    vg_regions = []
    for j in range(0, len(vg_genes)):
        vg_regions.append(string.upper(vg_genes[j].seq))

    vd_regions = []
    for j in range(0, len(vd_genes)):
        vd_regions.append(string.upper(vd_genes[j].seq))
        
    ja_regions = []
    for j in range(0, len(ja_genes)):
        ja_regions.append(string.upper(ja_genes[j].seq))

    jb_regions = []
    for j in range(0, len(jb_genes)):
        jb_regions.append(string.upper(jb_genes[j].seq))

    jg_regions = []
    for j in range(0, len(jg_genes)):
        jg_regions.append(string.upper(jg_genes[j].seq))

    jd_regions = []
    for j in range(0, len(jd_genes)):
        jd_regions.append(string.upper(jd_genes[j].seq))

    return va_regions, vb_regions, vg_regions, vd_regions, ja_regions, jb_regions, jg_regions, jd_regions

def build_keyword_tries(seqs):

    builder = AcoraBuilder()
    for i in range(0,len(seqs)):
        builder.add(str(seqs[i])) # Add all V tags to keyword trie

    key = builder.build()
    return key

def engine(rc, recid,
           seq_count, assigned_count, Nseqs,
           vb_key, jb_key,
           vb_seqs, jb_seqs,
           half1_vb_seqs, half2_vb_seqs,
           half1_jb_seqs, half2_jb_seqs,
           jump_to_end_vb, jump_to_start_jb,
           vb_regions, jb_regions,
           half1_vb_key, half2_vb_key,
           half1_jb_key, half2_jb_key,
           out_file,
           error0_count, error1_count,
           barcode_seq, barcode_qual
           ):

    ### Open .txt file created at the start of analysis
    stemplate = Template('$v $j $del_v $del_j $nt_insert $seqid $barcode $barqual') # Creates stemplate, a holder, for f. Each line will have the 5 variables separated by a space
    found_seq_match = 0
    found_v_match = 0
    found_j_match = 0
    
    print len(jb_seqs)
    
    hold_v = vb_key.findall(str(rc))
    hold_j = jb_key.findall(str(rc))
    
    v_match, temp_end_v, found_v_match,  error0_count, error1_count  = v_analysis( str(rc), hold_v, vb_seqs, half1_vb_seqs, half2_vb_seqs, jump_to_end_vb, vb_regions, half1_vb_key, half2_vb_key, error0_count, error1_count )
    j_match, temp_start_j, found_j_match = j_analysis( str(rc), hold_j, jb_seqs, half1_jb_seqs, half2_jb_seqs, jump_to_start_jb, jb_regions, half1_jb_key, half2_jb_key )
    
    if v_match != None and j_match != None:
      if get_v_deletions( str(rc), v_match, temp_end_v, vb_regions ) \
	and get_j_deletions( str(rc), j_match, temp_start_j, jb_regions ) \
	  and found_v_match == 1 \
	    and found_j_match == 1 :
	      [end_v, deletions_v] = get_v_deletions( str(rc), v_match, temp_end_v, vb_regions )
	      [start_j, deletions_j] = get_j_deletions( str(rc), j_match, temp_start_j, jb_regions )
	      f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(rc[end_v+1:start_j])+str(','), seqid = str(recid)+str(','), barcode = barcode_seq+str(','), barqual = barcode_qual )
	      if not ((temp_end_v - jump_to_end_vb[v_match]) + deletions_v) > ((temp_start_j -  deletions_j) + jump_to_start_jb[j_match]) or \
		not deletions_v > (jump_to_end_vb[v_match] - len(vb_seqs[v_match])) or \
		  not deletions_j > jump_to_start_jb[j_match]:
		  print >> out_file, f_seq # Write to out_file (text file) the classification of the sequence
		  assigned_count += 1
		  found_seq_match = 1

    return assigned_count, seq_count, Nseqs, found_seq_match, error0_count, error1_count 

def v_analysis( rc, hold_v, v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v, v_regions, half1_v_key, half2_v_key, error0_count, error1_count ):

    # rc is a string of record.seq as input

    v_match = None

    if hold_v:                
        v_match = v_seqs.index(hold_v[0][0]) # Assigns V
        temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
        if get_v_deletions( rc, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
            [ end_v, deletions_v] = get_v_deletions( rc, v_match, temp_end_v, v_regions )
        found_v_match = 1
        error0_count += 1
    else:
        found_v_match = 0
        hold_v1 = half1_v_key.findall(rc)
        hold_v2 = half2_v_key.findall(rc)
        for i in range(len(hold_v1)):
            indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
            for k in indices:
                if len(v_seqs[k]) == len(rc[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                    if lev.hamming( v_seqs[k], str(rc)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                        v_match = k
                        temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                        found_v_match += 1
                        error1_count += 1
        for i in range(len(hold_v2)):
            indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
            for k in indices:
                if len(v_seqs[k]) == len(rc[hold_v2[i][1]-v_half_split:hold_v2[i][1]-v_half_split+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                    if lev.hamming( v_seqs[k], rc[hold_v2[i][1]-v_half_split:hold_v2[i][1]+len(v_seqs[k])-v_half_split] ) <= 1:
                        v_match = k
                        temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1 # Finds where the end of a full V would be
                        found_v_match += 1
                        error1_count += 1

    if v_match is not None:
        return v_match, temp_end_v, found_v_match, error0_count, error1_count
    else:
        return [None, None, None, error0_count, error1_count]

def j_analysis( rc, hold_j, j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j, j_regions, half1_j_key, half2_j_key ):

    j_match = None
    
    if hold_j:
        j_match = j_seqs.index(hold_j[0][0]) # Assigns J
        temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
        if get_j_deletions( rc, j_match, temp_start_j, j_regions ):
            [ start_j, deletions_j] = get_j_deletions( rc, j_match, temp_start_j, j_regions )
        found_j_match = 1
    else:
        found_j_match = 0
        hold_j1 = half1_j_key.findall(rc)
        hold_j2 = half2_j_key.findall(rc)
        
        for i in range(len(hold_j1)):
            indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
            for k in indices:
                if len(j_seqs[k]) == len(rc[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                    if lev.hamming( j_seqs[k], rc[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                        j_match = k
                        temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                        found_j_match += 1
        for i in range(len(hold_j2)):
            indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
            for k in indices:
                if len(j_seqs[k]) == len(rc[hold_j2[i][1]-j_half_split:hold_j2[i][1]-j_half_split+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                    if lev.hamming( j_seqs[k], rc[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                        j_match = k
                        temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                        found_j_match += 1

    if j_match is not None:
        return j_match, temp_start_j, found_j_match
    else:
        return [None, None, None]

def get_v_deletions( rc, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence rc
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in rc.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= function_temp_end_v < len(rc):
        if str(v_regions_cut[v_match])[pos] == str(rc)[function_temp_end_v] and str(v_regions_cut[v_match])[pos-1] == str(rc)[function_temp_end_v-1] and str(v_regions_cut[v_match])[pos-2] == str(rc)[function_temp_end_v-2]:
            is_v_match = 1
            deletions_v = -pos - 1
            end_v = function_temp_end_v
        else:
            pos -= 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        return []

def get_j_deletions( rc, j_match, temp_start_j, j_regions_cut ):
    # This function determines the number of J deletions in sequence rc
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in rc.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j < len(str(rc)):
        if str(j_regions_cut[j_match])[pos] == str(rc)[function_temp_start_j] and str(j_regions_cut[j_match])[pos+1] == str(rc)[function_temp_start_j+1] and str(j_regions_cut[j_match])[pos+2] == str(rc)[function_temp_start_j+2]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        return []

def get_v_tags(file_v, half_split):
    import string
    
    v_seqs = [] # Holds all V tags
    jump_to_end_v = [] # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip("\n") # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(string.split(elements)[0]) # Adds elements in first column iteratively
        jump_to_end_v.append(int(string.split(elements)[1])) # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    import string
    
    j_seqs = [] # Holds all J tags
    jump_to_start_j = [] # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(string.split(elements)[0])
        jump_to_start_j.append(int(string.split(elements)[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]

def get_distinct_clones( handle, handle_results, with_count=False ):

    ## LOOKS THROUGH TEXT FILE OF CLASSIFIERS AND WRITES NEW FILE CONTAINING ALL DISTINCT CLASSIFIERS, OPTIONALLY WITH COUNT OF ALL DISTINCT CLASSIFIERS
    ## with_count=True writes file with counts of all distinct classifiers
    
    from string import Template
    import collections as coll
    from operator import itemgetter, attrgetter

    write_to = open(str(handle_results)+'.txt', "w")

    if with_count == True:
        stemplate = Template('$count $element')
        d = coll.defaultdict(int)
        for line in handle:
            classifier = line.rstrip("\n")
            elements = classifier.split(",")
            del elements[5:len(elements)]
            elements = str(elements)[1:-1]
            if elements in d:
                d[elements] += 1
            else:
                d[elements] = 1
        d_sorted = sorted(d.items(), key=itemgetter(1), reverse=True)
        for k in d_sorted:
            kcount = k[1]
            z = k[0].split(',')
            details = z[0].strip("'")+', '+z[1].strip("' ")+', '+z[2].strip("' ")+', '+z[3].strip("' ")+', '+z[4].strip("' ")
            f_seq = stemplate.substitute( count = str(kcount)+str(','), element = details )
            print >> write_to, f_seq
    else:
        stemplate = Template('$element')
        d = coll.defaultdict(int)
        for line in handle:
            classifier = line.rstrip("\n")
            elements = classifier.split(",")
            del elements[5:len(elements)]
            elements = str(elements)[1:-1]
            if elements in d:
                d[elements] += 1
            else:
                d[elements] = 1
        d_sorted = sorted(d.items(), key=itemgetter(1), reverse=True)
        for k in d_sorted:
            kcount = k[1]
            z = k[0].split(',')
            details = z[0].strip("'")+', '+z[1].strip("' ")+', '+z[2].strip("' ")+', '+z[3].strip("' ")+', '+z[4].strip("' ")
            f_seq = stemplate.substitute( element = details )
            print >> write_to, f_seq
    
    handle.close()
    write_to.close()

def get_translated_sequences( handle, handle_results, chain, species, with_outframe=False, fullaaseq=False ):

    ## TRANSLATES CLASSIFIERS TO AA SEQUENCES VIA THEIR NT SEQUENCE
    ## Default settings are -
    ## chain = "beta" or chain = "alpha"
    ## with_outframe=True or False: writes all aa seqeunces to file, including those that are out-of-frame (with stop codon symbol *)
    ## fullaaseq=True or False: True writes the whole V(D)J aa sequence to file, False, writes only the CDR3 region.

    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.Alphabet import generic_dna
    import string
    import re

    if species == 'human':
        handle_vb=open("human_TRBV_region.fasta","rU")
        handle_jb=open("human_TRBJ_region.fasta","rU")
        handle_va=open("human_TRAV_region.fasta","rU")
        handle_ja=open("human_TRAJ_region.fasta","rU")
        handle_vg=open("human_TRGV_region.fasta","rU")
        handle_jg=open("human_TRGJ_region.fasta","rU")
        handle_vd=open("human_TRDV_region.fasta","rU")
        handle_jd=open("human_TRDJ_region.fasta","rU")
    elif species == 'mouse':
        handle_vb=open("mouse_TRBV_region.fasta","rU")
        handle_jb=open("mouse_TRBJ_region.fasta","rU")
        handle_va=open("mouse_TRAV_region.fasta","rU")
        handle_ja=open("mouse_TRAJ_region.fasta","rU")
        handle_vg=open("mouse_TRGV_region.fasta","rU")
        handle_jg=open("mouse_TRGJ_region.fasta","rU")
        handle_vd=open("mouse_TRDV_region.fasta","rU")
        handle_jd=open("mouse_TRDJ_region.fasta","rU")
    
    vb_raw = list(SeqIO.parse(handle_vb, "fasta"))
    handle_vb.close()
    jb_raw = list(SeqIO.parse(handle_jb, "fasta"))
    handle_jb.close()
    va_raw = list(SeqIO.parse(handle_va, "fasta"))
    handle_va.close()
    ja_raw = list(SeqIO.parse(handle_ja, "fasta"))
    handle_ja.close()
    vg_raw = list(SeqIO.parse(handle_vg, "fasta"))
    handle_vg.close()
    jg_raw = list(SeqIO.parse(handle_jg, "fasta"))
    handle_jg.close()
    vd_raw = list(SeqIO.parse(handle_vd, "fasta"))
    handle_vd.close()
    jd_raw = list(SeqIO.parse(handle_jd, "fasta"))
    handle_jd.close()

    vb_regions = []
    for i in range(0,len(vb_raw)):
        vb_regions.append(string.upper(vb_raw[i].seq))

    jb_regions = []
    for i in range(0,len(jb_raw)):
        jb_regions.append(string.upper(jb_raw[i].seq))

    va_regions = []
    for i in range(0,len(va_raw)):
        va_regions.append(string.upper(va_raw[i].seq))

    ja_regions = []
    for i in range(0,len(ja_raw)):
        ja_regions.append(string.upper(ja_raw[i].seq))

    vg_regions = []
    for i in range(0,len(vg_raw)):
        vg_regions.append(string.upper(vg_raw[i].seq))

    jg_regions = []
    for i in range(0,len(jg_raw)):
        jg_regions.append(string.upper(jg_raw[i].seq))

    vd_regions = []
    for i in range(0,len(vd_raw)):
        vd_regions.append(string.upper(vd_raw[i].seq))

    jd_regions = []
    for i in range(0,len(jd_raw)):
        jd_regions.append(string.upper(jd_raw[i].seq))

    write_to = open( str(handle_results)+'.txt', "w")

    if chain == "alpha":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 8:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ',''))
            elif len(classifier) == 7:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''
                
            if delv != 0:
                used_v = va_regions[v][:-delv]
            elif delv == 0:
                used_v = va_regions[v]

            if delj != 0:
                used_j = ja_regions[j][delj:]
            elif delj == 0:
                used_j = ja_regions[j]

            seq = str(used_v + ins + used_j)
            start = (len(seq)-1)%3
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G',str(aaseq)) and re.findall('C',str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper+4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3

    if chain == "beta":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 8:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ',''))
            elif len(classifier) == 7:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''

            if delv != 0:
                used_v = vb_regions[v][:-delv]
            elif delv == 0:
                used_v = vb_regions[v]

            if delj != 0:
                used_j = jb_regions[j][delj:]
            elif delj == 0:
                used_j = jb_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq)%3+2
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G',str(aaseq)) and re.findall('C',str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper+4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3

    if chain == "gamma":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 8:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ',''))
            elif len(classifier) == 7:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''

            if delv != 0:
                used_v = vg_regions[v][:-delv]
            elif delv == 0:
                used_v = vg_regions[v]

            if delj != 0:
                used_j = jg_regions[j][delj:]
            elif delj == 0:
                used_j = jg_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq)%3+2
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G',str(aaseq)) and re.findall('C',str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper+4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3

    if chain == "delta":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 8:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ',''))
            elif len(classifier) == 7:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''

            if delv != 0:
                used_v = vd_regions[v][:-delv]
            elif delv == 0:
                used_v = vd_regions[v]

            if delj != 0:
                used_j = jd_regions[j][delj:]
            elif delj == 0:
                used_j = jd_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq)%3+2
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G',str(aaseq)) and re.findall('C',str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper+4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3
            
    handle.close()
    write_to.close()