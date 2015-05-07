# FindDs.py v1.0

# Mines decombinator insert data for the presence of Diversity regions (either TRBD or TRDD)
  # Outputs results weighted by frequency of reads
    # e.g if there were 2 unique rearrangements, one with TRBD1 present once and one with TRBD2 present twice
      # this scenario would be 33.33% TRBD1, 66.67% TRBD2

from __future__ import division
import collections as coll
import difflib as dl
import sys
#import os
import matplotlib.pyplot as plt
import numpy as np
import time 

fontsize = 12
plt.rcParams.update({'font.size': fontsize})

savepath = "/home/jme/TCR/WRITE_UP/THESIS/WorkingPlots/" + time.strftime("%Y %m %d") + " "

len_threshold = 7

if len(sys.argv) < 3:
  print "Please supply a filename containing vDCR-decombined data and chain (b/d), e.g. python CollapseTCRs.py FILENAME.freq b"
  print "You may also supply an optional extra field specifying threshold length of D region allowed for assignation"
  print "e.g. python CollapseTCRs.py FILENAME.freq b 8 (default is 7)"
  sys.exit()
else:
  dcrfilename = str(sys.argv[1])
  dcrfile = open(dcrfilename, "rU")
  chain = sys.argv[2]
  if len(sys.argv) == 4:
    len_threshold = int(sys.argv[3])
  if len(sys.argv) > 4:
    print "Please supply a filename containing vDCR-decombined data and chain (b/d), e.g. python CollapseTCRs.py FILENAME.freq b"
    print "You may also supply an optional extra field specifying threshold length of D region allowed for assignation"
    print "e.g. python CollapseTCRs.py FILENAME.freq b 8 (default is 7)"
    sys.exit()    
    
    
    
  
if chain not in ['b', 'd']:
  print "The only valid chains are \'b\' and \'d'\, as only TRB and TRD have D genes"

bdnames = ['TRBD1','TRBD2','TRBD2*2'] 

bd = ['GGGACAGGGGGC', 'GGGACTAGCGGGGGG', 'GGGACTAGCGGGAGGG']

ddnames = ['TRDD1', 'TRDD2', 'TRDD3']
dd = ['GAAATAGT', 'CCTTCCTAC', 'ACTGGGGGATACG']


uniq_dcr_count = 0
tot_dcr_count = 0 
uniq_with_d = 0
uniq_with_dd = 0
tot_with_d = 0
tot_with_dd = 0

# dictionaries to store the frequency of all the associated V genes 
  # also initialise with empty values for all possible genes (i.e. number of extended tags)
bd1vs = coll.Counter()
bd1js = coll.Counter()
bd2vs = coll.Counter()
bd2js = coll.Counter()

# dictionaries to store length distribution (only for NON TANDEM genes)
unidlens = coll.Counter()	# those we can't ID
d1lens = coll.Counter()
d2lens = coll.Counter()

for i in range(63):
  bd1vs[str(i)] = 0
  bd2vs[str(i)] = 0

for i in range(13):
  bd1js[str(i)] = 0
  bd2js[str(i)] = 0
  
unid = 0 	# counter for number of unidentifiable D genes
canid = 0	# and those we can ID

for dcr in dcrfile:
  # loop through all records in Decombined file
  
  bits = dcr.rstrip().split(", ")
  insert = bits[4]

  uniq_dcr_count += 1
  tot_dcr_count += int(bits[5])
  
  possds = coll.Counter()
  
  found = False 	# check that only count each matching dcr once
  double_chk = False # need to check haven't already counted a DCR (otherwise tandems get counted twice)
  
  for g in range(len(vars()[chain + "d"])):
    # look for all known D genes of that chain
    s = dl.SequenceMatcher(None, insert, vars()[chain + "d"][g])
    match = s.find_longest_match(0, len(insert), 0, len(vars()[chain + "d"][g]))
    
    if match[2] >= len_threshold:
      possds[vars()[chain + "dnames"][g]] = match[2]
      
      #print insert, vars()[chain + "d"][g], vars()[chain + "dnames"][g], insert[match[0]:match[0]+match[2]], bits[0]
    
    if possds and found==False: 
      uniq_with_d += 1
      tot_with_d += int(bits[5])
      found = True
      
    #sys.exit()
	
    # check for tandem Ds, i.e. VDDJ rearrangements
    
    if chain == "b" and double_chk == False:
      if ('TRBD1' in possds.keys() and 'TRBD2' in possds.keys()) or ('TRBD1' in possds.keys() and 'TRBD2*2' in possds.keys()): 
	uniq_with_dd += 1
	tot_with_dd += int(bits[5])
	double_chk = True
	#print dcr.rstrip()
	#sys.exit()
    elif chain == "d" and double_chk == False:
      if len(possds.keys()) > 1:
	uniq_with_dd += 1
	tot_with_dd += int(bits[5])
	double_chk = True
    #if len(possds.keys()) >1: sys.exit()
      ##sys.exit()
      
  # count lengths of Ds    
  if len(possds.keys()) == 1 and double_chk == False:
    if possds.most_common()[0][0] == 'TRBD1':
      d1lens[possds.most_common()[0][1]] += int(bits[5])
    elif possds.most_common()[0][0] in ['TRBD2', 'TRBD2*2']:
      d2lens[possds.most_common()[0][1]] += int(bits[5])
  elif len(possds.keys()) > 1:
    if sorted([possds.most_common()[0][0], possds.most_common()[1][0]]) == ['TRBD2', 'TRBD2*2']:
      # many times TRDB2 will get two equal hits
      d2lens[possds.most_common()[0][1]] += int(bits[5])     
    elif possds.most_common()[0][1] == possds.most_common()[1][1]:      
      unidlens[possds.most_common()[0][1]] += int(bits[5])
    else:
      if possds.most_common()[0][0] == 'TRBD1':
	d1lens[possds.most_common()[0][1]] += int(bits[5])
      elif possds.most_common()[0][0] in ['TRBD2', 'TRBD2*2']:
	d2lens[possds.most_common()[0][1]] += int(bits[5])     

      
  # want to record which Js different Ds associate with (only for singles)
  
  if found == True and double_chk == False:
       
    if len(possds.keys()) == 1:
      canid += 1
      tempgene = possds.most_common()[0][0]
    
    else:
      firstmatchlen = possds.most_common()[0][1]
      secondmatchlen = possds.most_common()[1][1]
     
      if firstmatchlen == secondmatchlen:
	unid += 1
	continue

      else:
	canid += 1
	tempgene = possds.most_common()[0][0]
      
    if tempgene == 'TRBD1':
      bd1vs[bits[0]] += int(bits[5])
      bd1js[bits[1]] += int(bits[5])
    elif tempgene == 'TRBD2' or tempgene == 'TRBD2*2':
      bd2vs[bits[0]] += int(bits[5])
      bd2js[bits[1]] += int(bits[5])
    
pc_uniq_with_d = uniq_with_d/uniq_dcr_count * 100
pc_uniq_with_dd = uniq_with_dd/uniq_dcr_count * 100

pc_tot_with_d = tot_with_d/tot_dcr_count * 100
pc_tot_with_dd = tot_with_dd/tot_dcr_count * 100

#print '{0:,}'.format(uniq_dcr_count), "records analysed, looking for TR" + chain.upper() + "D regions with a threshold of", str(len_threshold), "nucleotides"
#print '{0:,}'.format(uniq_with_d), "contained detectable D region =", str(round(pc_uniq_with_d, 2)), "%"
#print '{0:,}'.format(uniq_with_dd), "contained detectable tandem regions =", str(round(pc_uniq_with_dd, 2)), "%"

print '{0:,}'.format(tot_dcr_count), "records analysed, looking for TR" + chain.upper() + "D regions with a threshold of", str(len_threshold), "nucleotides"
print '{0:,}'.format(tot_with_d), "contained detectable D region =", str(round(pc_tot_with_d, 2)), "%"
print '{0:,}'.format(tot_with_dd), "contained detectable tandem regions =", str(round(pc_tot_with_dd, 2)), "%"





#print '{0:,}'.format(), ""

# how I ran this on all files 
  # I first gzipped all the CD4/8 sorted data as well as the v1 bleeds for those as have them
  # therefore get one sample analysed per donor, just for whole repertoires
  # then put all into a table
#for x in {5..13}; do echo $x; for i in *beta*HV*q; do echo $i; python FindDs.py $i b $x >> $x.trbd; done; done 

#also apply to alphas, to get a feel for false positive rate
#jme@jaybuntu:/media/jme/SAMSUNG/ThesisAnalysis/BLR$ for x in {5..13}; do echo $x; for i in *alpha*HV*q; do echo $i; python FindDs.py $i b $x >> $x.atrbd; done; done 


sys.exit()

###### plotting D length distributions

# ensure all dictionaries have right number of keys

for x in range(len_threshold, max(d2lens.keys())+1):
  unidlens[x] += 0
  d1lens[x] += 0
  d2lens[x] += 0

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

ax.bar(d2lens.keys(), d2lens.values(), bottom=0, color="black")
ax.bar(d1lens.keys(), d1lens.values(), bottom=d2lens.values(), color="white")
ax.bar(unidlens.keys(), unidlens.values(), bottom=[x+y for x,y in zip(d2lens.values(), d1lens.values())], color="grey")

plt.show()




######### plotting gene usage

# TRBJ pairing

pts1 = []
pts2 = []

for i in range(13):
  pts1.append(bd1js[str(i)]/(sum(bd1js.values())+sum(bd2js.values())))
  pts2.append(bd2js[str(i)]/(sum(bd1js.values())+sum(bd2js.values())))


xs = np.arange(13)
width=.4
trbj = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7']

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

ax.set_xticks([x+.45 for x in range(13)])
ax.set_xticklabels(trbj, rotation=90)

plt.bar(xs, pts1, width, color="black", label="TRBD1")
plt.bar(xs+width, pts2, width, color="lightgray", label="TRBD2")

plt.xlim(0,13)
plt.legend(loc="upper left", prop={'size':fontsize})
plt.ylabel("Proportion of rearrangements")


plt.savefig(savepath + str(len_threshold) + "_" + dcrfilename.split('_')[2].split('.')[0] + "_TRBDJ_Pairing.svg", bbox_inches='tight')
#plt.show()
plt.close()

###############

# TRBV pairing

pts1 = []
pts2 = []

trbv = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 57, 58]
trbj = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
trbvnam = ['TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2', 'TRBV20-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1', 'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5', 'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-8', 'TRBV6-9', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV7-9', 'TRBV9', 'TRBV17', 'TRBV23-1', 'TRBV5-3', 'TRBV5-7', 'TRBV6-7', 'TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV21-1', 'TRBV22-1', 'TRBV5-2', 'TRBV7-5']

for i in trbv:
  pts1.append(bd1vs[str(i)]/(sum(bd1vs.values())+sum(bd2vs.values())))
  pts2.append(bd2vs[str(i)]/(sum(bd1vs.values())+sum(bd2vs.values())))

xs = np.arange(len(trbv))
width=.4


fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)

ax.set_xticks([x+.45 for x in range(len(trbv))])
ax.set_xticklabels(trbvnam, rotation=90)

plt.bar(xs, pts1, width, color="black", label="TRBD1")
plt.bar(xs+width, pts2, width, color="lightgray", label="TRBD2")
plt.ylabel("Proportion of rearrangements")
plt.xlim(0,len(trbv))
plt.legend(loc="upper left", prop={'size':fontsize})

plt.savefig(savepath + str(len_threshold) + "_" + dcrfilename.split('_')[2].split('.')[0] + "_TRBDV_Pairing.svg", bbox_inches='tight')
#plt.show()
plt.close()















