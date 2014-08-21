## Executes functions from functionsv2.x
## according to arguments parsed below

import sys, argparse, os
import DecombinatorFunctionsV2_2 as f
import Plotting as p
import ShortReadDecombinator as ShortReads
import warnings

warnings.filterwarnings("ignore") ## Supress warning from Biopython stated when translating of sequences where n mod 3 != 0.

parser = argparse.ArgumentParser(description='Decombinator v2.2')
parser.add_argument('-i','--input', help='Enter the path to the input fastq file you wish to analyse', required=True)
parser.add_argument('-o','--output', help='Enter the name you wish to call the output results files', required=True)
parser.add_argument('-rev','--reversecomplement', help='Enter True or False for whether a search on reverse complement of sequences is also desired', required=False, default=True, type=bool)
parser.add_argument('-b','--barcode', help='Enter True or False for whether sequences contain barcodes', required=False, default=False, type=bool)
parser.add_argument('-bs1','--barcodebegin1', help='Enter integer defining the start of the first barcode region', required=False, type=int )
parser.add_argument('-be1','--barcodefinish1', help='Enter integer defining the end of the first barcode region', required=False, type=int )
parser.add_argument('-bs2','--barcodebegin2', help='Enter integer defining the start of the second barcode region', required=False, type=int )
parser.add_argument('-be2','--barcodefinish2', help='Enter integer defining the end of the second barcode region', required=False,type=int )
parser.add_argument('-p','--withplots', help='Enter True or False', required=False, default=False, type=bool)
parser.add_argument('-sh','--shortreads', help='Enter True or False', required=False, default=False, type=bool)
parser.add_argument('-c','--count', help='Enter True or False', required=False, default=False, type=bool)
parser.add_argument('-of','--outofframe', help='Enter True or False', required=False, default=False, type=bool)
parser.add_argument('-f','--fullsequence', help='Enter True or False', required=False, default=False, type=bool)
parser.add_argument('-ch','--chaintype', help='Enter alpha, beta, gamma, delta or all', required=False, default='all')
parser.add_argument('-s','--speciestype', help='Enter human or mouse', required=False, default='human')
args = vars(parser.parse_args())

inputfile = args['input']
outputfile = args['output']
revsearch = args['reversecomplement']
barcoding = args['barcode']
barcodestart1 = args['barcodebegin1']
barcodeend1 = args['barcodefinish1']
barcodestart2 = args['barcodebegin2']
barcodeend2 = args['barcodefinish2']
include_plots = args['withplots']
forshortreads = args['shortreads']
withcount = args['count']
outframe = args['outofframe']
fullseq = args['fullsequence']
chain = args['chaintype']
species = args['speciestype']

newpath = f.create_folder(outputfile)

if forshortreads == True:
    vfile = os.getcwd()+'/human_TRBV_region.fasta'
    jfile = os.getcwd()+'/human_TRBJ_region.fasta'
    v_key, j_key, v_regions, j_regions = ShortReads.setup(vfile,jfile)
    #infile = str(inputfile)
    fileid = str(outputfile)
    param_set = [10, 2, 1400, 1.05]
    ShortReads.analyse_file( inputfile, newpath, fileid, v_key, j_key, v_regions, j_regions, param_set)
else:
    f.analysis( inputfile, outputfile, with_reverse_complement_search=revsearch, barcode=barcoding, barcodestart1=barcodestart1, barcodeend1=barcodeend1, barcodestart2=barcodestart2, barcodeend2=barcodeend2, newpath=newpath, omitN=True, chain=chain, species=species)

chains = ['alpha','beta','delta','gamma']
for chain in chains:
    if os.stat(newpath+outputfile+'_'+chain+'.txt').st_size != 0: # if the file is non-empty, i.e. if TcRchain seqs were found...
        print 'Getting distinct clones for TcR'+chain
        f.get_distinct_clones( open(newpath+outputfile+'_'+chain+'.txt', "rU"), handle_results=newpath+str('distinct_clones')+'_'+chain,with_count=withcount )
        print 'Translating sequences for TcR'+chain
        f.get_translated_sequences( open(newpath+outputfile+'_'+chain+'.txt', "rU"), handle_results=newpath+str('translated_sequences')+'_'+chain, chain=str(chain), species=species, with_outframe=outframe, fullaaseq=fullseq )
        
if include_plots==True:
    print 'Plotting the results of the analysis'
    if forshortreads == True:
        if os.stat(newpath+outputfile+'_beta'+'.txt').st_size != 0: # if the file is non-empty, i.e. if TcRchain seqs were found...
            p.plot_v_usage( open(newpath+outputfile+'_beta'+'.txt', "rU"), chain='beta', species=species, savefilename = newpath+'Vusage')
            p.plot_j_usage( open(newpath+outputfile+'_beta'+'.txt', "rU"), chain='beta', species=species, savefilename = newpath+'Jusage')
            p.plot_del_v( open(newpath+outputfile+'_beta'+'.txt', "rU"), savefilename = newpath+'Vdels')
            p.plot_del_j( open(newpath+outputfile+'_beta'+'.txt', "rU"), savefilename = newpath+'Jdels')
            p.plot_vj_joint_dist( open(newpath+outputfile+'_beta'+'.txt', "rU"), chain=str(chain), species=species, savefilename = newpath+'VJusage')
            p.plot_insert_lengths( open(newpath+outputfile+'_beta'+'.txt', "rU"), savefilename = newpath+'InsertLengths')
    else:    
        chains = ['alpha','beta','delta','gamma']
        for chain in chains:
            if os.stat(newpath+outputfile+'_'+chain+'.txt').st_size != 0: # if the file is non-empty, i.e. if TcRchain seqs were found...
                p.plot_v_usage( open(newpath+outputfile+'_'+chain+'.txt', "rU"), chain=chain, species=species, savefilename = newpath+'Vusage'+chain )
                p.plot_j_usage( open(newpath+outputfile+'_'+chain+'.txt', "rU"), chain=chain, species=species, savefilename = newpath+'Jusage'+chain )
                p.plot_del_v( open(newpath+outputfile+'_'+chain+'.txt', "rU"), savefilename = newpath+'Vdels'+chain)
                p.plot_del_j( open(newpath+outputfile+'_'+chain+'.txt', "rU"), savefilename = newpath+'Jdels'+chain)
                p.plot_vj_joint_dist( open(newpath+outputfile+'_'+chain+'.txt', "rU"), chain=str(chain), species=species, savefilename = newpath+'VJusage'+chain)
                p.plot_insert_lengths( open(newpath+outputfile+'_'+chain+'.txt', "rU"), savefilename = newpath+'InsertLengths'+chain)
