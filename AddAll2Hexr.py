### HISTORY/PURPOSE ###
# Developed as an alternative to AddN12to6.py
# Adds all R2 hexamers to R1 file (and actually produces a third file) BEFORE decombining
# Loops through two fastqs simultaneously, and outputs a third = R1 (which has our VDJ) with the first 6 bases of R2 (R2 hex) stuck at the the very front

### RUNNING/INPUT ###
# Runs off command line input of 3 file names, where first 2 are (Illumina encoded) fastq files, and third is an output fastq
# Run: $ python AddAllR2Hex.py read1.fastq read2.fastq output.fq
# Run on all fastq pairs in directory: $ 
# for i in *R1*; do j=${i/R1/R2}; k=$(echo $i | cut -d "_" -f1); echo $k; python AddAllR2Hex.py $i $j $k.fq; done
# Expect ~4-5 minutes per million reads on an average spec machine (seems to do about 1/4 of a million reads per min)
# INPUT FASTQ FILES MUST BE PAIRED END! Script will work on any 2 fastqs, but will only produce a meaningful result in paired-end files of equal length and corresponding read positions



### OUTPUT ###
# A fastq file (of specified name) that (should!) contain the VDJ recombination, where the first 12 bases are the random barcode

from Bio import SeqIO
from time import time, clock
from itertools import izip
import sys

filename = ""

if (len(sys.argv) <> 4):
  print "Please supply 2 input and one output file names (i.e. python AddAllR2Hex.py read1.fastq read2.fastq output.fq)"
  sys.exit()
else:
  fq1file = str(sys.argv[1])
  fq2file = str(sys.argv[2])
  outfq = str(sys.argv[3])

fq1 = SeqIO.parse(open(fq1file), "fastq")
fq2 = SeqIO.parse(open(fq2file), "fastq")

outfile = open(outfq, "w")

count = 0

t0 = time() # Begin timer

for record1, record2 in izip(fq1, fq2):
  
  ### For non-standard Illumina encoded fastqs, might need to change which fields are carried into fq_* vars
  
  fq_id = record1.id
  fq_seq = record2.format("fastq").split('\n')[1][:6] + str(record1.seq)
  fq_qual = record2.format("fastq").split('\n')[3][:6] + record1.format("fastq").split('\n')[3]
  
  new_record = str("@" + fq_id + "\n" + fq_seq + "\n+\n" + fq_qual + "\n")
  
  outfile.write(new_record)
  
  count += 1

outfile.close()

timed = time() - t0

print count, 'reads processed from', fq1file, 'and', fq2file, 'and output into', outfq
print '\t\t\t\t\t\t\t\t\tTook', round(timed,2), 'seconds'