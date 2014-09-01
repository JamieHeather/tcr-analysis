# CollapseTCRs.py v1.2
# James M. Heather, August 2014, UCL

##################
### BACKGROUND ###
##################

# Makes use of random barcode sequences to error-correct high-throughput sequencing TCR repertoire data
# Random nucleotides (N) added at known positions to amplicons prior to amplification
  # Therefore these can act as molecular barcodes of original cDNA molecules, allowing us to cancel out PCR amplification/error
# Consists of two parts:
  # First, error correction: 
    # The script collapses TCR sequences within a barcode, correcting for sequencing and PCR error
  # Next, PCR-amplification mitigation:
    # The barcodes that contribute to a particular TCR are clustered and counted, giving a more accurate frequency value 
    
# Also note that in the lab we have created a neater, more rationally laid out version of this program (via Katharine Best)
  # This will hopefully feature in a later publication 
  # However this is still the version I used to generate the data for my thesis and the HIV-TCR paper

##################
###### INPUT #####
##################

# Takes the .n12 output from vDCR.py (which contains TCR rearrangement and unique barcode information)
# Run: python CollapseTCRs.py FILENAME.n12 

##################
##### OUTPUT #####  
##################
    
# Outputs a .freq file = standard 5-part TCR identifier of collapsed TCRs, plus a sixth, collapsed-frequency field

##################
### DEFINITIONS ##  
##################

# dcr: TCR 5-part identifier as defined by any version of Decombinator, or sometimes Decombinator itself (DeCombinatoR, see?)
# dcr_etc: lines of an n12 file that pass the barcode filter then get stored in this '|' delimited string for processing
# .n12: File produced from verbose Decombinator (vDCR), which includes TCR and barcode sequence and quality information
# fdf: frequency dcr file, i.e. the output file into which the decombinator identifiers plus frequency get written
# clash: different TCRs (DCRs) associating with same barcode
# dodex/dodec/dodecamer: all refer to the dodecameric random nucleotide barcode sequence, introduced before PCR

##################
#### PACKAGES ####  
##################

from __future__ import division
import collections as coll
import sys
import Levenshtein as lev
import re
from operator import itemgetter
from time import time, clock
import json
import signal
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

if (len(sys.argv) <> 2):
  print "Please supply a filename containing vDCR-decombined data, e.g. python CollapseTCRs.py FILENAME.n12"
  sys.exit()
else:
  seqdcrfilename = str(sys.argv[1])
  seqdcrfile = open(seqdcrfilename, "rU")
  
start_time = time() # set time at start 

###################################################################
################# DEFINE FUNCTIONS, COUNTERS ETC ##################
###################################################################  

def breakdown(etc):                                
  # Used to break a given dcr_etc (i.e. what is stored in a given line of a .n12 file) into its components
  # Splits on '|', to avoid breaking up within dcr identifiers or fastq quality strings
  
  return re.findall(r"[\w,\<\=\>\-\;\:\?\/\.\@\# ]+", str(etc))
  
  # breakdown[0] = dcr / [1] = seq / [2] = qual / [3] = id
  
def fq(etc):
  # Rebuilds a dcr_etc instance back into a full fastq record
  
  fq_id = breakdown(etc)[3]
  fq_str = breakdown(etc)[1]
  fq_qual = breakdown(etc)[2]
  fastq_string = "@%s\n%s\n+\n%s\n" % (fq_id, fq_str, fq_qual)
  record = SeqIO.read(StringIO(fastq_string), "fastq")
  return (record)

global ldist    # FIX - do I need this?
def ldist(one, two):
  # Return the Levenshtein distance of two strings
  return lev.distance(one, two)

def cluster():
  # Acts on DCR-collapsed dictionary, i.e. that which has already been error-corrected
  # Looks within a given TCR (DCR) at how many barcode clusters there are, which is used to estimate frequency of original cDNA
    # A cluster is all those barcodes that are within a given number of edits from all other members
    # For each TCR, it takes all the associated barcodes, arranged decreasing by order, and looks at each in turn
    # The first (biggest) barcode gets put in a cluster 
    # The next is compared to that; if it's <= threshold (default set to 3) it gets clustered with that sequence
    # If >= threshold, it starts a new cluster. Continue until all barcodes are assigned to a cluster
  # Therefore within a cluster all members are within a certain distance of all other members of that group
  
  # Note that while the threshold seems quite high, the regions of the reads that we employed for our barcode
    # (in this iteration of the protocol) were the first 6 nucleotides of reads 1 and 2
    # Unfortunately both of these sections suffer from poor quality data, for different reasons
    
  # Note that this function was produced with the very generous help of Katharine Best
  
  distance_threshold = 3        # The number of edits allowed to be included in the same cluster
  
  freq_dcr_file = open(fdf, "w")
  freq_dcr_file.close()
  
  global col_count
  global count_barcodes_kept
  
  with open(fdf, "a") as a:
    for x in dcr_collapsed:
      
      # If there is only one associated barcode, it must have a frequency of 1
      if len(dcr_collapsed[x]) == 1:
        count_barcodes_kept += 1
        result = str(x) + ", 1\n"
        
      else:
        tracking = coll.Counter()
        clustered = coll.Counter()
        
        # Use dictionaries to keep track of which barcodes have been seen, or tracking new barcodes as needed
        for k,v in dcr_collapsed[x].most_common():      
          
          if tracking[k] == 0:
            ref_bc = k
            tracking[k] = 1
            clustered[ref_bc] += v
            
            for k2,v2 in dcr_collapsed[x].most_common():
              
              if tracking[k2] == 0:
                if ldist(ref_bc, k2) <= distance_threshold:
                  tracking[k2] = 1
                  clustered[ref_bc] += v2
                
        count_barcodes_kept += sum(clustered.values())
        result = str(x) + ", " + str(len(clustered)) + "\n"
        
      a.write(result)
      col_count += 1      
  
### Counters

count_clash = 0                 # clash = seemingly different TCRs (or DCRs at least) having same dodecamer
count_dodec = 0                 # count number dodex    
count_dcrseq = 0                # count number of sequences

count_match = 0
count_same_len_mismatch = 0
count_same_len_discarded = 0
count_diff_len_mismatch = 0
count_diff_len_discarded = 0

count_biggest_fail = 0          # count if biggest clone not properly assigned
count_long = 0

count_barcode_fail = 0          # check for barcodes containing low quality score calls

global count_barcodes_kept

count_barcodes_kept = 0         # number of different barcodes that actually end up contributing


### Dictionaries

def dodex():
  return coll.Counter()

global dcr_collapsed

dcr_collapsed = coll.defaultdict(dodex)

clone = coll.defaultdict(list)                  # dd to hold dodecamer {dcr/seq}

count_dcrs = coll.Counter()                     # counter to hold unique dcrs, so we know what fraction survive the collapse

clustered = coll.defaultdict(int)               # holds eventual clustered TCR-DCR counts


###################################################################
################### GENERATE CLONE DICT ###########################
###################################################################


  ### COMMA DELIMITED POSITIONS IN THIS DATA ARE, IN ORDER: ###
  
  # V -[0]- J -[1]- Vdel -[2]- Jdel -[3]- insert -[4]- ID -[5]- TCRseq -[6]- TCRqual -[7]- barcode -[8]- barcode qual -[9]- 
  

print "\nReading", str(seqdcrfilename), "into dictionary..."

t0 = time() # Begin timer

count_input_lines = 0

for line in seqdcrfile:
  
  count_input_lines += 1

  comma = [m.start() for m in re.finditer(',', line)]           #define commas, from which we define all our 

  dcr = str(line[:comma[4]])
  
  ID = str(line[comma[4]+2:comma[5]])
  
  seq = str(line[comma[5]+2:comma[6]])
  
  qual = str(line[comma[6]+2:comma[7]])
  
  barcode = str(line[comma[7]+2:comma[8]])
  
  # Want to only allow assignations if their barcode is of sufficiently high quality
    # This is the most important string to be sure of (given its importance, length and propensity to be poor quality in R2)
  
  barcode_qual = str(line[comma[8]+2:])         
  
  barcode_etc = "barcode" + "|" + barcode + "|" + barcode_qual + "|" + ID 
  
  barcode_fq = fq(barcode_etc)
  
  barcode_check = [s for s in barcode_fq.letter_annotations.values()[0] if s < 20]
  
  barcode_sum_qual = sum(barcode_fq.letter_annotations.values()[0])
  
  #######################
  ###  BARCODE CHECK! ### 
  #######################
  
  # Only allows barcodes sequences through if two requirements met: 
    # they contain a maximum of 1 base with a quality below 20
    # the sum total of their qualities exceeds 360 (equivalent average base call being >= Q30) 
  
  if len(barcode_check) > 1 or barcode_sum_qual < 360:  
    count_barcode_fail += 1

  else:
    dcr_etc = dcr + "|" + seq + "|" + qual + "|" + ID 
    clone[barcode].append(dcr_etc)
    count_dcrs[dcr] += 1
    
timed = time() - t0

print '\t\t\t\t\t\t\tTook', round(timed,2), 'seconds'


###################################################################
#################### FIND BIGGEST CLONE ###########################
###################################################################

print "Collapsing clones, removing error..."

t0 = time() # Begin timer

  
for d in clone:                         # loop through each different dodecamer
 
  #print "--------------------------------------------------------------------------"
  count_dodec += 1
 
  # Loops through all clones twice, once to establish the most common, and then again to collapse with reference to that
  
  dcr_count = coll.Counter()
  
  for c in clone[d]:            # first loop, find most common clone
  
    count_dcrseq += 1           # count each dcrseq, ie each dcr assignation
    
    clone_seq = breakdown(c)[0] + "|" + breakdown(c)[1]
    
    dcr_count[clone_seq] += 1 
  
  biggest_max_size = max(dcr_count.values())            # size of biggest clone (slightly better way!)
  
  biggest_clones_list = [k for k in dcr_count.keys() if dcr_count[k] == biggest_max_size]               
    # produce a list of clones that have the same frequency as the biggest

  if len(biggest_clones_list) > 1:      
    # when there is more than 1 DCR that has the same frequency as the biggest, i.e. 2+ in joint 1st place
    
    biggest_avg_qual = 0
    biggest_dcrseq = ""
    
    for cc in clone[d]:                 # loop to find sequences of biggest clones            
      
      for b in biggest_clones_list:     # then look to see...
        
        cc_gaps = [m.start() for m in re.finditer("\|", cc)]
        
        if b == cc[:cc_gaps[1]]:        
	  # ... if this is the instance that got put in biggest_clones_list (based on its DCR-str combo)
          
          record = fq(cc)               # turn it back into a fastq
          
          avg_qual = sum(record.letter_annotations.values()[0]) / len(record.letter_annotations.values()[0])    
	    # use biopython to interpret the quality scores, of the new record, get an avg
           
          if avg_qual > biggest_avg_qual:
            biggest_avg_qual = avg_qual
            biggest_dcrseq = b                  # by the end of the loop, this will be the biggest
                         
          # NB if there are two clones that both the same frequency and the same overall average quality
	    # then the script will just keep whichever was checked last in the loop
	  # However the odds of this happening are very small, thus this is very unlikely to be a problem with large clones

      if biggest_avg_qual == 0:
        count_biggest_fail += 1
        
  else:                 # or if there's an obvious top clone
    biggest_dcrseq = biggest_clones_list[0]



###################################################################
############### COMPARE ALL AGAINST BIGGEST #######################
###################################################################

  # Having found the largest clone within a barcode, we then assume that to be the 'genuine' clone for that barcode
  # Then all other clones within the barcode are compared against it (using the nucleotide sequence)
    # If they are within a given threshold then they are assumed to be errors from that genuine clone
    # If they are above that threshold they are assumed to be either too error prone to include, or genuine barcode clashes
  # NB Effectively this just throws out all clones that aren't the biggest
    # However early testing showed that almost all reads do fall within the threshold; very few are thrown out
    # Theoretically secondary rounds of collapsing could be included for those clones that are thrown out
      # However given our throughput this was merited to really not be worth it; we just don't get enough clashes
      
  for x in clone[d]:                    # loop to decide whether to collapse/keep each new clone
    
    clone_seq = breakdown(x)[0] + "|" + breakdown(x)[1]
    
    dcr = breakdown(x)[0]
    seq = breakdown(x)[1]

    proto_dcr = breakdown(biggest_dcrseq)[0]    
    proto_seq = breakdown(biggest_dcrseq)[1]
    
  ### Perform various checks, and compare the current clone to the prototypical (or 'biggest') clone
  
    if len(seq) > 130:                  # sequence length sanity check
      count_long += 1
      
      ### RECORD DISMISSED - NOT CARRIED OVER ###
       
    elif seq == proto_seq:           # exact match with prototypical
      
      count_match += 1
      
      ##### OUTPUT THIS RECORD! ######
      
      dcr_collapsed[proto_dcr][d] += 1
      
    elif len(seq) == len(proto_seq):
      
      count_same_len_mismatch += 1
      count_clash += 1
      
      dist = lev.distance(seq, proto_seq)  # calculate levenshtein distance between the biggest sequence and the current
      
      pc_dist = (dist/len(seq)) * 100           # calc percentage difference
      
      if pc_dist > 20:  
        
        count_same_len_discarded += 1
        
        ### RECORD DISMISSED - NOT CARRIED OVER ###
        
      else:
        
       count_match += 1
      
        ##### OUTPUT THIS RECORD! ######
      
       dcr_collapsed[proto_dcr][d] += 1
                    
    elif len(proto_seq) <> len(seq):             # if sequences are different lengths
    
      count_diff_len_mismatch += 1
      count_clash += 1
      
      dist = lev.distance(seq, proto_seq)  # calculate levenshtein distance between the biggest sequence and the current
      
      pc_dist = (dist/len(seq)) * 100         
      
      if pc_dist <= 20:
        
        count_match += 1

        #### OUTPUT THIS RECORD! ######
      
        dcr_collapsed[proto_dcr][d] += 1
      
      elif pc_dist > 20:
        
        # need to check if the distance is being artificially inflated by DCR finding a tag up/down-stream of the genuine ones
        # can try trimming either end of the longer sequence to see if the distance then drops low enough
        
        dist1 = 0
        dist2 = 0
        dist3 = 0
        dist4 = 0
        
        if len(seq) > len(proto_seq):
          
          diff = len(seq) - len(proto_seq)
          
          dist1 = lev.hamming(seq[:len(proto_seq)], proto_seq)   # trim 3'
          
          dist2 = lev.hamming(seq[diff:], proto_seq)             # trim 5'
          
          if dist1/len(proto_seq) * 100 <= 20 or dist2/len(proto_seq) * 100 <= 20:
                    
            count_match += 1

            #### OUTPUT THIS RECORD! ######
            
            dcr_collapsed[proto_dcr][d] += 1
            
          else:

            count_diff_len_discarded += 1
            
            ### RECORD DISMISSED - NOT CARRIED OVER ###            
                   
        elif len(proto_seq) > len(seq):
          
          diff = len(proto_seq) - len(seq)
           
          dist3 = lev.hamming(proto_seq[:len(seq)], seq)         # trim 3'
          
          dist4 = lev.hamming(proto_seq[diff:], seq)             # trim 5' 
          
          if dist3/len(proto_seq) * 100 <= 20 or dist4/len(proto_seq) * 100 <= 20:
                    
            count_match += 1

            #### OUTPUT THIS RECORD! ######
            
            dcr_collapsed[proto_dcr][d] += 1
            
          else:

            count_diff_len_discarded += 1
            
            ### RECORD DISMISSED - NOT CARRIED OVER ###            

timed = time() - t0

print '\t\t\t\t\t\t\tTook', round(timed,2), 'seconds'


###################################################################
############### CLUSTER BARCODES INTO FREQUENCIES #################
###################################################################

# Given the propensity of R2 barcodes to be poor quality, we have relatively stringent criterita for barcode inclusion
  # Can't afford to be as lenient as we were with the actual TCR sequences
# Only include barcodes if they occur more than once (filtering out a great many)
  # unless of course there is only one barcode for a given dcr_collapsed
# Then add barcodes to clusters where each member is within 3 edits of any other
  # The number of these clustered barcodes then represents a more conservative estimate of cDNA frequency than number of reads

print 'Clustering barcodes into frequencies...'

global fdf

fdf = seqdcrfilename.split(".")[0]+str(".freq")

global col_count
col_count = 1

time_elapsed = time() - start_time

t0 = time() # Begin timer    

try:
  cluster()
except Exception, msg:
    print str(msg)

timed = time() - t0

print '\t\t\t\t\t\t\tTook', round(timed,2), 'seconds'

print '{0:,}'.format(len(dcr_collapsed)), "unique clones output to", fdf

##################
## OUTPUT STATS ##
##################

# Here are a number of optional statistics (primarily relating to error-correction) that will print to stdout if desired

##### Comment this out if you want to see the statistics
sys.exit()
#####

print "\n--------------------------------------------------------------\n"

print '{0:,}'.format(count_dcrseq), "sequences processed:", '{0:,}'.format(len(count_dcrs)), "different assigned TCRs associated with", '{0:,}'.format(count_dodec), "dodecamers\n"

print '{0:,}'.format(count_match), "clones match prototypical exactly \t\t=", round(((count_match/count_dcrseq)*100), 2), "% of total sequences"

print '{0:,}'.format(count_clash), "clashes detected pre-processing \t\t\t=", round(((count_clash/count_dcrseq)*100), 2), "% of total sequences\n"

print count_same_len_mismatch, "clones were same length, but not 100% identical \t=", round(((count_same_len_mismatch/count_clash)*100), 2), "% of total clashes"
print "Of these,", count_same_len_discarded, "were discarded as they were more than 20% different to the prototypical clone (genuine dodec clashes?)\n"

print count_diff_len_mismatch, "clones were not the same length as prototypical \t=", round(((count_diff_len_mismatch/count_clash)*100), 2), "% of total clashes"
print "Of these,", count_diff_len_discarded, "were discarded as they were more than 20% different to the prototypical clone (genuine dodec clashes?)\n"

print count_biggest_fail, "clones failed to assign a prototypical clone"
print count_long, "clones failed length sanity check (> 130 nt)\n"
print '{0:,}'.format(count_barcode_fail), "sequences discarded due to low quality barcode sequences"

print '{0:,}'.format(count_barcodes_kept), "different barcodes contributed to the final frequency counts"

sys.exit()

