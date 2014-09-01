# AddAllR2Hex.py v1.2
# James M. Heather, January 2014, UCL

##################
### BACKGROUND ###
##################

# Used to turn paired-end fastq data from our barcoded TCR-amplification protocol into correct format for vDCR.py
# Adds all R2 random hexamers to start of R1 file, meaning those reads start with a random 12-mer barcode
# Loops through two fastqs simultaneously, and outputs a third, in format:
  # |-R2N6-|-R1N6-|-----------R1-V(D)J--------------------|
    # Where R1N6 = read 1 random hexamer, R2N6 = read 2 random hexamer, and R1-V(D)J is the TCR rearrangement found in R1

##################
###### INPUT #####
##################

# Runs off command line input of 3 file names, where first 2 are (Illumina encoded) fastq files, and third is an output fastq:
  # python AddAllR2Hex.py read1.fastq read2.fastq output.fq
# Run on all fastq pairs in directory: 
  # for i in *R1*; do j=${i/R1/R2}; k=$(echo $i | cut -d "_" -f1); echo $k; python AddAllR2Hex.py $i $j $k.fq; done
# Expect ~4-5 minutes per million reads on an average spec machine (seems to do about 1/4 of a million reads per min)

# NB: Input files must be paired-end (and correctly ordered)
  # Script will work on any 2 fastqs, but will only produce a meaningful result on properly paired end data
    # (that has been produced using our protocol, or a closely analagous one)

##################
##### OUTPUT #####  
##################

# A fastq file (of specified name) that contains the VDJ recombination, with the 12-mer random barcode at the start of the read

##################
#### PACKAGES ####  
##################

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
  
  ### For non-standard Illumina encoded fastqs, one might need to change which fields are carried into fq_* vars
  
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