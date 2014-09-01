# CDR3ulator.py v1.9
# James M. Heather, August 2014, UCL

##################
### BACKGROUND ###
##################

# Take any decombined data and output the functional CDR3s only
  # CDR3s must be: in-frame; lacking-stop codons, and run from a conserved cysteine to FGXG motif (or appropriate alternatives)
# Originally built on PlotDCR.py (itself modified from dcr v1.4), then algorithm swapped
  # Now uses CDR3 detection based on Katharine's functions.py script
  # Uses defined conserved cysteine residues, and allows for atypical TRAJ CDR3-ending motifs
  
# Makes use of functions originally developed in Katharine Best's functions.py script and Niclas Thomas' Decombinator (v1.2)
# Note that this version provides the capabilities to generate CDR3s from all genes covered in the extended Decombinator tags
  # Many of these genes will never generate CDR3s from this code regardless of whether they're included in CDR3s
  # This is typically due to V genes that lack the conserved C residue, or contain stop codons upstream of it.
  # These have been left in, both to provide a single location that describes the whole heterogeneity of the prototypical alpha/beta genes

##################
###### INPUT #####
##################

# Takes any text file in comma-space (", ") delimited decombinator format
  # 5-part TCR identifier must come first, followed by any additional fields, such as frequency
# TCR chain must either be specified explicitly (with "a" or "b", or contained within the input filename
# Run: python CDR3ulator.py beta_FILE.txt
# Or: python CDR3ulator.py FILE.freq a
# Note providing an explicit chain character overrides any chain string present in filename

##################
##### OUTPUT #####  
##################

# Users have options of two output, comma-delimited formats:
  # '.cdr3' files, which consist of the unique productively-rearranged CDR3s from the original file, with frequencies
  # '.dcrcdr3' files, which contains the five-part Decombinator index before the CDR3 it encodes and its frequency
  # The choice of which file format is used is decided by altering the 'dcr_output' variable
  # Note that output files will likely be shorter than the input file, due to multiple TCR sequences encoding the same CDR3
# Users also have the option to include the final three residues of the CDR3 in the output, or to stop at the phenylalanine 

# The script also outputs some simple statistics about the number of productive rearrangements, and the frequency of different functionality genes
# NB: See IMGT for further information on definitions of functionality and productivity:
  # http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html

##################
#### PACKAGES ####  
##################

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import string
import re
from Bio import SeqIO
import sys
import re
import collections as coll
import os
import urllib2

##################
## USER OPTIONS ##
##################

dcr_output = False	# False => .cdr3 files / True => .dcrcdr3 files
include_gxg = False	
        
###################
#### FUNCTIONS ####
###################

def findfile(testfile):

  try:
      testopen = open(str(filename),"rU")
      testopen.close()
  except:
      print 'Cannot find the file you specified. Please try again'
      sys.exit()

def import_gene_information(species, chain):            

  # Runs first: reads in V and J gene sequence and name data (from fasta files)
    # and positions of conserved cysteine residues in V genes (from separate files)
  
  # If files cannot be found in local directory, script looks for them online at GitHub
  
  # NB that a number of psuedogenes have no officially designated conserved C (or indeed a 3' C at all)
    # Where possible, the nearest suitable C residue is used, where not an arbitrary position of 0 is given
    # Somewhat moot, as most psuedogenes contain a number of stop codons and thus cannot produce productive rearrangements 
    
  if os.path.isfile("exthuman_TR" +string.upper(chain[0])+ "V_region.fasta"):
    vfile = "exthuman_TR" +string.upper(chain[0])+ "V_region.fasta"
    with open(vfile, 'r') as f:
      vfasta = list(SeqIO.parse(f, 'fasta'))
  else:
    onlineV = "https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TR" + string.upper(chain[0]) + "V_region.fasta"

    vfasta = list(SeqIO.parse(urllib2.urlopen(onlineV), 'fasta'))
      
  vregions = [str(string.upper(item.seq)) for item in vfasta]  
  vnames = [str(string.upper(item.id).split("|")[1]) for item in vfasta]    

  if os.path.isfile("exthuman_TR" +string.upper(chain[0])+ "J_region.fasta"):
    vfile = "exthuman_TR" +string.upper(chain[0])+ "J_region.fasta"
    with open(vfile, 'r') as f:
      jfasta = list(SeqIO.parse(f, 'fasta'))
  else:
    onlineJ = "https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/exthuman_TR" + string.upper(chain[0]) + "J_region.fasta"
    jfasta = list(SeqIO.parse(urllib2.urlopen(onlineJ), 'fasta'))

  jregions = [str(string.upper(item.seq)) for item in jfasta]  
  jnames = [str(string.upper(item.id).split("|")[1]) for item in jfasta]    

  if os.path.isfile("TR" +string.upper(chain[0])+ "V_ConservedC.txt"):
    with open("TR" +string.upper(chain[0])+ "V_ConservedC.txt", 'r') as f:
      vconservedc = [int(line.rstrip('\n')) for line in f]
  else:
    onlineVC = "https://raw.githubusercontent.com/JamieHeather/tcr-analysis/master/TR" + string.upper(chain[0]) + "V_ConservedC.txt"
    vconservedc = [int(line.rstrip('\n')) for line in urllib2.urlopen(onlineVC)]
    
  return vregions, jregions, vnames, jnames, vconservedc

def get_cdr3(dcr, chain, vregions, jregions, vconservedc):

  # Checks the productivity of a given DCR-assigned rearrangement 
  # Returns a 1 if productive, 0 if not 
  # NB: A productively rearranged receptor does not necessarily mean that it is the working receptor used in a cell
    # It could be a silenced chain that isn't used, or could have inactivating mutations upstream of the sequenced region
  
  # 0.5 Set up check variables
  
  # Boolean productivity checks that CDR3s must pass 
  in_frame = 0
  no_stop = 0
  found_c = 0
  found_fgxg = 0
  
  # CDR3-defining positions 
  start_cdr3 = 0
  end_cdr3 = 0 
  
  # 1. Rebuild whole nucleotide sequence from Decombinator assignment
  classifier_elements = dcr.split(', ')
  v = int(classifier_elements[0])
  j = int(classifier_elements[1])
  vdel = int(classifier_elements[2])
  jdel = int(classifier_elements[3])
  ins_nt = classifier_elements[4]

  if vdel == 0:
    v_used = vregions[v]
  else:
    v_used = vregions[v][:-vdel]
  j_used = jregions[j][jdel:]

  nt = ''.join([v_used, ins_nt, j_used])
  
  # 2. Check whether whole rearrangement is in frame
  if (len(nt)-1) % 3 == 0:
    in_frame = 1
  else:
    return "Out of frame"

  # 3. Translate
  aa = str(Seq(nt, generic_dna).translate())

  # 4. Check for stop codons
  if '*' not in aa:
    no_stop = 1
  else:
    return "Contains stop codon"

  # 5. Check for conserved cysteine in the V gene
  if aa[vconservedc[v]-1] == 'C':
    found_c = 1
    start_cdr3 = vconservedc[v]-1
  elif chain == "beta" and v in [45, 56]: # Allow for TRBV17 and TRBV26, which use Y instead of C to start CDR3s
    if aa[vconservedc[v]-1] == 'Y':
      found_c = 1
      start_cdr3 = vconservedc[v]-1
    
  else:
    return "No conserved cysteine"
  
  # 5.5 Having found conserved cysteine, only need look downstream to find other end of CDR3
  downstream_c = aa[start_cdr3:]

  # 6. Check for presence of FGXG motif (or equivalent)
  if chain == 'alpha':      

    # TRAJs get a little more interesting with their conserved sequences than the other genes
      # All TRAJ FGXGs are at the -11 position apart from one
        # TRAJ59 (DCR 58) is truncated in its 3', and thus its FGXG is at -9:-5
      # All TRAJS use FGXG as their CDR3 ending motif, apart from the following:
        #TRAJ16 (DCR 6) = FARG
        #TRAJ33 (DCR 22) = WGAG
        #TRAJ38 (DCR 26) = WGLG
        #TRAJ35 (DCR 54) = CGSG
        #TRAJ51 (DCR 55) = FGKE
        #TRAJ55 (DCR 56) = WGKG
        #TRAJ61 (DCR 60) = FGAN      
    
    if j <> 58:
      site = downstream_c[-11:-7]
      
      if re.findall('FG.G', site):
	if include_gxg == True:	
	  end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 4 
	else:
	  end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1
      
      else: # Allow for the non-canonical FGXGs in TRAJ      

	awkward_ajs = [6, 22, 26, 54, 55, 56, 60]
	alt_aj_motifs = ['FARG', 'WGAG', 'WGLG', 'CGSG', 'FGKE', 'WGKG', 'FGAN']
	
	if j in awkward_ajs and site == alt_aj_motifs[awkward_ajs.index(j)]:
	  if include_gxg == True:	
	    end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 4       
	  else:
	    end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1
   
	else:
	  return "No conserved FGXG"
	  
    else: # TRAJ59
      site = downstream_c[-9:-5]     
      
      if re.findall('FG.G', site):
	if include_gxg == True:	
	  end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 4 
	else:
	  end_cdr3 = len(downstream_c) - 11 + start_cdr3 + 1 
    
  elif chain == 'beta':
    # All F TRBJ FGXGs being at -10 position
      # TRBJ2-2P is at -8:-4 (although I am unaware of any evidence in our own or others' data that it gets recombined)
    
    site = downstream_c[-10:-6]
    
    if re.findall('FG.G', site):
      if include_gxg == True:
	end_cdr3 = len(downstream_c) - 10 + start_cdr3 + 4 
      else:
        end_cdr3 = len(downstream_c) - 10 + start_cdr3 + 1

    elif j == 13: # TRBJ2-2Ps
      site = downstream_c[-8:-4]
      
      if re.findall('LGGG', site):
	if include_gxg == True:
	  end_cdr3 = len(downstream_c) - 8 + start_cdr3 + 4 
	else:
	  end_cdr3 = len(downstream_c) - 8 + start_cdr3 + 1      

  else:
    return "No conserved FGXG"
  
  return aa[start_cdr3:end_cdr3]    

####################################
######## CHECK INPUT FILES #########
####################################

if (len(sys.argv) == 3):

  # If chain is explicitly provided:
  
  filename = str(sys.argv[1])
  findfile(filename)
  inputchain = str(sys.argv[2])
  if inputchain == "a" or inputchain == "alpha":
    chain = "alpha"
  elif inputchain == "b" or inputchain == "beta":
    chain = "beta"
  else:
    print "Can't assign chain. Indicate chain with either \'a\' or \'b\', or include in filename"
    sys.exit
 
elif (len(sys.argv) == 2):
  
  # If chain read from input file name:
  filename = str(sys.argv[1])
  a = "alpha"
  b = "beta"
  if a.upper() in filename.upper():
    chain = "alpha"
  elif b.upper() in filename.upper():
    chain = "beta"
  else:
    print "Can't assign chain. Make sure decombined file name contains either \'alpha\' or \'beta\'"
    sys.exit()
 
else:
  
  print "Incorrect input command. Please supply either a filename containing a chain, or both a file and a chain letter"
  print "e.g.: python CDR3ulator.py FILE.txt a"
  print "or: python CDR3ulator.py FILE_beta.freq"       

####################################
########## EXTRACT CDR3s ###########
####################################

info = import_gene_information("human", chain)
  # Returns: vregions, jregions, vnames, jnames, vconservedc

infile = open(filename, "rU")

line_count = 0
F_count = 0

# Count non-productive rearrangments
fail_count = coll.Counter()
fails = ["Out of frame", "Contains stop codon", "No conserved cysteine", "No conserved FGXG"]

# Store and count CDR3s
if dcr_output == False:
  cdr3_count = coll.Counter()
elif dcr_output == True:
  dcr_cdr3_count = coll.Counter()

# Count the number of F, ORF and P genes used for both productive and non-productive rearrangements
pVf = 0
pVorf = 0
pVp = 0
pJf = 0
pJorf = 0
pJp = 0

npVf = 0
npVorf = 0
npVp = 0
npJf = 0
npJorf = 0
npJp = 0

for line in infile:
  
  line_count += 1
  
  comma = [m.start() for m in re.finditer(',', line)]           

  in_dcr = str(line[:comma[4]])
  
  cdr3 = get_cdr3(in_dcr, chain, info[0], info[1], info[4])
  
  frequency = int(line[comma[4]+2:].rstrip())
  
  dcr_cdr3 = in_dcr + ":" + cdr3
  
  v = int(line[:comma[0]])
  j = int(line[comma[0]+2:comma[1]])
  
  if cdr3 not in fails:
  
    F_count += 1    
    
    if dcr_output == False:
      cdr3_count[cdr3] += frequency
    elif dcr_output == True:
      dcr_cdr3_count[dcr_cdr3] += frequency
	
    # Count the number of number of each type of gene functionality (by IMGT definitions)
  
    # TRAV: 0 -> 46 = F; 47 = ORF; 48 -> 55 = P
    # TRBV: 0 -> 44 = F; 45 -> 50 = ORF; 51 -> 62 = P
	# NB 45 + 56 USE Y INSTEAD OF C
    # TRAJ: 0 -> 49 = F; [50,51,52,53,54,57,58,60] = ORF; [55, 56, 59] = P
    # TRBJ: 0 -> 12 = F; 13 = ORF
	
    # NB Many P V genes lack conserved cysteine residues: if there is a neighbouring C, that has been used
      # Similarly many P J genes lack conserved FGXG motifs: nearest conserved motif used
    
    if chain == "alpha":
      if v <= 46:
	pVf += 1
      elif v == 47:
	pVorf += 1
      elif v >= 48:
	pVp += 1    
      
      if j <= 49:
	pJf += 1
      elif j in [50,51,52,53,54,57,58,60]:
	pJorf += 1
      elif j in [55, 56, 59]:
	pJp += 1
      
    elif chain == "beta":
      if v <= 44:
	pVf += 1
      elif v >= 45 and v <= 50:
	pVorf += 1
      elif v >= 51:
	pVp += 1    
	
      if j <= 12:
	pJf += 1
      elif j == 13:    
	pJorf += 1
    
  else:
    
    fail_count[cdr3] += 1
    
    # And again for the non-productively rearranged receptors

    if chain == "alpha":
      if v <= 46:
	npVf += 1
      elif v == 47:
	npVorf += 1
      elif v >= 48:
	npVp += 1    
      
      if j <= 49:
	npJf += 1
      elif j in [50,51,52,53,54,57,58,60]:
	npJorf += 1
      elif j in [55, 56, 59]:
	npJp += 1
      
    elif chain == "beta":
      if v <= 44:
	npVf += 1
      elif v >= 45 and v <= 50:
	npVorf += 1
      elif v >= 51:
	npVp += 1    
	
      if j <= 12:
	npJf += 1
      elif j == 13:    
	npJorf += 1    
	
 # Uncomment to output information on each line read in (for debugging purposes)    
  #print v, info[2][v], j, info[3][j], "\t", cdr3 
  
if dcr_output == True:

  outfilename = filename.split(".")[0]+".dcrcdr3"
  outfile = open(outfilename, "w")
    
  for x in dcr_cdr3_count:
    
    outtext = x + ", " + str(dcr_cdr3_count[x])
    
    print >> outfile, outtext
      
  infile.close()
  outfile.close()
  
elif dcr_output == False:

  outfilename = filename.split(".")[0]+".cdr3"
  outfile = open(outfilename, "w")

  for x in cdr3_count:
    
    outtext = x + ", " + str(cdr3_count[x])
    
    print >> outfile, outtext
      
  infile.close()
  outfile.close()

NF_count = sum(fail_count.values())
        
###################
##### RESULTS #####
###################

print "Reading", str(line_count), "Decombinator-assigned rearrangements from", str(filename) + ", and writing out to", str(outfilename), "\n"

print '{0:,}'.format(F_count), "productive rearrangements detected" 
print "\tV gene usage:\t" + '{0:,}'.format(pVf), "F;\t" + '{0:,}'.format(pVorf), "ORF;\t" +  '{0:,}'.format(pVp), "P"
print "\tJ gene usage:\t" + '{0:,}'.format(pJf), "F;\t" + '{0:,}'.format(pJorf), "ORF;\t" + '{0:,}'.format(pJp), "P\n"

print '{0:,}'.format(NF_count), "non-productive rearrangements detected"
print "\tV gene usage:\t" + '{0:,}'.format(npVf), "F;\t" + '{0:,}'.format(npVorf), "ORF;\t" +  '{0:,}'.format(npVp), "P"
print "\tJ gene usage:\t" + '{0:,}'.format(npJf), "F;\t" + '{0:,}'.format(npJorf), "ORF;\t" + '{0:,}'.format(npJp), "P\n"

if (F_count + NF_count) <> line_count:
  print "\nError detected: Sum of productive and non-productive sorted sequences does not equal total number of input sequences"
  