# PlotDCR.py v1.3
# James M. Heather, August 2014, UCL

##################
### BACKGROUND ###
##################

# Take any Decombinator file and output graphs of V/J usage, insertions and deletions
# Allow visual comparison of differently processed DCR files without having to re-run the analysis
# Core code is taken straight from Decombinator v1.4 (see https://github.com/uclinfectionimmunity/Decombinator/)
  # Separating analysis and plotting code allows processing of decombined data before plotting
# This version also provides extra functionality in that it allows plotting of the additional genes covered by the extended tags
  # i.e. this version will allow plotting of the prototypical versions of ORF and P genes, as well as the F
    # See vDCR.py and CDR3ulator.py scripts for more information (https://github.com/JamieHeather/tcr-analysis)

##################
###### INPUT #####
##################

# Should take any text file in comma-space (", ") delimited Decombinator format
  # Assumes 5-part identifier makes up the first 5 fields, and any frequency resides in sixth field
# Also note that for this implementation of the pipeline only works on human alpha/beta TCRs
# Chain must be specified in filename or given explicitly
# e.g. run 'python PlotDCR.py filename_beta.freq'
  # or 'python PlotDCR.py filename.txt a'
  
# Note that users can specify which tag set (original or extended) they wish to use for plotting
  # Set the 'default_tags' variable to "original" or "extended"
  # Original tags overrule extended, i.e. one can run original tags on Decombinator results from extended
    # Users may wish to do so as a large number of the extended tag genes will never appear without changing upstream code

##################
##### OUTPUT #####  
##################

# A new folder containing all plots and graphs: V and J gene usage and deletions, insert lengths and V/J pairing 
# Output files are named after the original input file

import sys
import os
import platform

# Here users can record which set of Decombinator tags have been used to create the data
default_tags = "original"

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
       

def plot_v_usage( handle, chain, savefilename="Vusage", order="frequency"):

###########################################################################################################################
############     ###  ######   ##        #####  #######  ##################################################################
############  ##  ##  #####  #  ####  #########  #####  ###################################################################
############     ###  ####  ###  ###  ##########  ###  ####################################################################
############  ######  #####  #  ####  ###########  #  #####################################################################
############  ######     ###   #####  ############   ######################################################################
###########################################################################################################################
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    if chain=="alpha":
	if tag_set == "original":
	  tags = open("tags_trav.txt", "rU")
	elif tag_set == "extended":
	  tags = open("exttags_trav.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            if (tag_set == "original" and int(elements.split(',')[0]) >= 46) :
	      continue	              
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order=="frequency":        
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_v)
            percent_usage_v = [0]*num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            
	    if tag_set == "original":
	      gene_list_v = ('V1-1','V1-2','V10','V12-1','V12-2','V12-3','V13-1','V13-2','V14/D4','V16','V17','V18','V19','V2','V20','V21','V22','V23/D6','V24','V25','V26-1','V26-2','V27','V29/DV5','V3','V30','V34','V35','V36/DV7','V38-1','V38-2/DV8','V39','V4','V40','V41','V5','V6','V7','V8-1','V8-2/8-4','V8-3','V8-6','V9-1','V9-2','DV1','DV2','DV3')
            v_linked = [0]*len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0]*len(sorted_v)
            v_percents = [0]*len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.yticks( pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)
            
        elif order=="chromosome":
            total = sum(freq_vector_v)
            fv = [0]*num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
	    if tag_set == "original":                
	      gene_list_v = ('V1-1','V1-2','V10','V12-1','V12-2','V12-3','V13-1','V13-2','V14/D4','V16','V17','V18','V19','V2','V20','V21','V22','V23/D6','V24','V25','V26-1','V26-2','V27','V29/DV5','V3','V30','V34','V35','V36/DV7','V38-1','V38-2/DV8','V39','V4','V40','V41','V5','V6','V7','V8-1','V8-2/8-4','V8-3','V8-6','V9-1','V9-2','DV1','DV2','DV3')
	      chromosome_order = [0,1,13,24,32,35,36,37,38,42,2,3,39,40,6,4,7,8,43,5,41,9,10,11,12,14,15,16,17,44,18,19,20,22,23,25,21,26,27,28,29,30,31,33,34,45,46]
	    elif tag_set == "extended":
	      gene_list_v = ('V1-1', 'V1-2', 'V10', 'V12-1', 'V12-2', 'V12-3', 'V13-1', 'V13-2', 'V14/DV4', 'V16', 'V17', 'V18', 'V19', 'V2', 'V20', 'V21', 'V22', 'V23/DV6', 'V24', 'V25', 'V26-1', 'V26-2', 'V27', 'V29/DV5', 'V3', 'V30', 'V34', 'V35', 'V36/DV7', 'V38-1', 'V38-2/DV8', 'V39', 'V4', 'V40', 'V41', 'V5', 'V6', 'V7', 'V8-1', 'V8-2/8-4', 'V8-3', 'V8-6', 'V9-1', 'V9-2', 'DV1', 'DV2', 'DV3', 'V8-7', 'V8-5', 'V11', 'V15', 'V28', 'V31', 'V32', 'V33', 'V37')
	      chromosome_order = [0, 1, 13, 24, 32, 35, 36, 37, 38, 42, 2, 49, 3, 39, 40, 6, 4, 48, 7, 8, 43, 50, 5, 41, 9, 10, 11, 12, 14, 15, 16, 17, 44, 18, 19, 20, 47, 22, 51, 23, 25, 52, 53, 54, 21, 26, 27, 28, 55, 29, 30, 31, 33, 34, 45, 46]
     
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize = 8)
            ax.set_xticks(ind+1*width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":
	if tag_set == "original":      
	  tags = open("tags_trbv.txt", "rU")
	elif tag_set == "extended":      
	  tags = open("exttags_trbv.txt", "rU")	  
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            if tag_set == "original" and int(elements.split(',')[0]) >= 45:
		continue
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order=="frequency":        
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_v)
            percent_usage_v = [0]*num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)

	    if tag_set == "original":                
	      gene_list_v = ('V10-1','V10-2','V10-3','V11-1','V11-2','V11-3','V12-3/V12-4','V12-5','V13','V14','V15','V16','V18','V19','V2','V20-1','V24-1','V25-1','V27','V28','V29-1','V3-1','V30','V4-1','V4-2','V4-3','V5-1','V5-4','V5-5','V5-6','V5-8','V6-1','V6-4','V6-5','V6-6','V6-8','V6-9','V7-2','V7-3','V7-4','V7-6','V7-7','V7-8','V7-9','V9')
            v_linked = [0]*len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0]*len(sorted_v)
            v_percents = [0]*len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.yticks( pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)
            
        elif order=="chromosome":
            total = sum(freq_vector_v)
            fv = [0]*num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)

	    if tag_set == "original":                
	      gene_list_v = ('V10-1','V10-2','V10-3','V11-1','V11-2','V11-3','V12-3/V12-4','V12-5','V13','V14','V15','V16','V18','V19','V2','V20-1','V24-1','V25-1','V27','V28','V29-1','V3-1','V30','V4-1','V4-2','V4-3','V5-1','V5-4','V5-5','V5-6','V5-8','V6-1','V6-4','V6-5','V6-6','V6-8','V6-9','V7-2','V7-3','V7-4','V7-6','V7-7','V7-8','V7-9','V9')
	      chromosome_order = [ 14, 21, 23, 26, 31, 24, 25, 37, 32, 38, 44, 0, 3, 1, 4, 33, 39, 27, 34, 28, 40, 29, 35, 41, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22 ]
	    elif tag_set == "extended":    
	      gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27', 'V28', 'V29-1', 'V3-1', 'V30', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9', 'V17', 'V23-1', 'V5-3', 'V5-7', 'V6-7', 'V7-1', 'V1', 'V12-1', 'V12-2', 'V21-1', 'V22-1', 'V26', 'V5-2', 'V7-5', 'V8-1', 'V8-2', 'VA', 'VB')
	      chromosome_order = [51, 14, 21, 23, 26, 31, 50, 24, 25, 37, 59, 57, 32, 38, 60, 47, 44, 0, 3, 52, 1, 4, 53, 33, 39, 27, 34, 58, 28, 49, 40, 29, 35, 41, 48, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 45, 12, 13, 15, 54, 55, 46, 16, 17, 61, 56, 62, 18, 19, 20, 22]     
	      
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+1*width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags.close()

    
def plot_j_usage( handle, chain="beta", savefilename="Jusage", order="frequency"):

###########################################################################################################################
############     ###  ######   ##        #####          ###################################################################
############  ##  ##  #####  #  ####  ############  #######################################################################
############     ###  ####  ###  ###  #######  ###  #######################################################################
############  ######  #####  #  ####  ########  #  ########################################################################
############  ######     ###   #####  #########   #########################################################################
###########################################################################################################################
 
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    if chain=="alpha":
        
	if tag_set == "original":        
	  tags = open("tags_traj.txt", "rU")
	elif tag_set == "extended":        
	  tags = open("exttags_traj.txt", "rU")	  
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            if (tag_set == "original" and int(elements.split(',')[1]) >= 49) :
	      continue	              
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order=="frequency":
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_j)
            percent_usage_j = [0]*num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)

	    if tag_set == "original":                
	      gene_list_j = ('J10','J11','J12','J13','J14','J15','J16','J17','J18','J20','J21','J22','J23','J24','J26','J27','J28','J29','J3','J30','J31','J32','J33','J34','J36','J37','J38','J39','J4','J40','J41','J42','J43','J44','J45','J46','J47','J48','J49','J5','J50','J52','J53','J54','J56','J57','J6','J7','J8','J9')
            j_linked = [0]*len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0]*len(sorted_j)
            j_percents = [0]*len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.yticks( pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_j)
            fj = [0]*num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)

	    if tag_set == "original":   
	      gene_list_j = ('J10','J11','J12','J13','J14','J15','J16','J17','J18','J20','J21','J22','J23','J24','J26','J27','J28','J29','J3','J30','J31','J32','J33','J34','J36','J37','J38','J39','J4','J40','J41','J42','J43','J44','J45','J46','J47','J48','J49','J5','J50','J52','J53','J54','J56','J57','J6','J7','J8','J9')
	      chromosome_order = [45,44,43,42,41,40,38,37,36,35,34,33,32,31,30,29,27,26,25,24,23,22,21,20,19,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,49,48,47,46,39,28,18]
	    elif tag_set == "extended":
	      gene_list_j = ('J10', 'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J20', 'J21', 'J22', 'J23', 'J24', 'J26', 'J27', 'J28', 'J29', 'J3', 'J30', 'J31', 'J32', 'J33', 'J34', 'J36', 'J37', 'J38', 'J39', 'J4', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J5', 'J50', 'J52', 'J53', 'J54', 'J56', 'J57', 'J6', 'J7', 'J8', 'J9', 'J1', 'J2', 'J19', 'J25', 'J35', 'J51', 'J55', 'J58', 'J59', 'J60', 'J61')
	      chromosome_order = [60, 59, 58, 57, 45, 44, 56, 43, 42, 41, 55, 40, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 27, 26, 25, 24, 54, 23, 22, 21, 20, 19, 17, 16, 15, 14, 53, 13, 12, 11, 10, 9, 52, 8, 7, 6, 5, 4, 3, 2, 1, 0, 49, 48, 47, 46, 39, 28, 18, 51, 50]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":
        
	if tag_set == "original":        
	  tags = open("tags_trbj.txt", "rU")
	elif tag_set == "extended":        
	  tags = open("exttags_trbj.txt", "rU")	  
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            if tag_set == "original" and int(elements.split(',')[1]) >= 12:
	      continue            
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order=="frequency":
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_j)
            percent_usage_j = [0]*num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)

	    if tag_set == "original":                   
	      gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
            j_linked = [0]*len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0]*len(sorted_j)
            j_percents = [0]*len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.yticks( pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_j)
            fj = [0]*num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)

	    if tag_set == "original":                   
	      gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
	      chromosome_order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
	    if tag_set == "extended":                   
	      gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7', 'J2-2P')
	      chromosome_order = [0, 1, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11, 12]

            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags.close()

def plot_del_v( handle, savefilename="Vdels"):

###########################################################################################################################
############     ###  ######   ##        #####  #######  ###     ###      ##  #############################################
############  ##  ##  #####  #  ####  #########  #####  ####  ##  ##  ######  #############################################
############     ###  ####  ###  ###  ##########  ###  #####  ###  #    ####  #############################################
############  ######  #####  #  ####  ###########  #  ######  ##  ##  ######  #############################################
############  ######     ###   #####  ############   #######     ###      ##    ###########################################
###########################################################################################################################

    import numpy as np
    import matplotlib.pyplot as plt
    import string

    if tag_set == "original":
      deletions_v = [0]*50
    elif tag_set == "extended":
      deletions_v = [0]*65      
      
    for line in handle:
        elements = line.rstrip("\n")
        deletions_v[int(elements.split(',')[2])] += 1

    total = sum(deletions_v)
    for i in range(len(deletions_v)):
        deletions_v[i] = deletions_v[i] / float(total)

    if tag_set == "original":
      ind = np.arange(50)
    elif tag_set == "extended":
      ind = np.arange(65)     
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_v, width, color='yellow')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of V germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_del_j( handle, savefilename="Jdels"):

###########################################################################################################################
############     ###  ######   ##        #####          #####     ###      ##  ############################################
############  ##  ##  #####  #  ####  ############  #########  ##  ##  ######  ############################################
############     ###  ####  ###  ###  #######  ###  #########  ###  #    ####  ############################################
############  ######  #####  #  ####  ########  #  ##########  ##  ##  ######  ############################################
############  ######     ###   #####  #########   ###########     ###      ##     #########################################
###########################################################################################################################
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_j = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_j[int(elements.split(',')[3])] += 1

    total = sum(deletions_j)
    for i in range(len(deletions_j)):
        deletions_j[i] = deletions_j[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_j, width, color='red')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of J germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_vj_joint_dist( handle, chain="beta", savefilename="VJusage" ):


###########################################################################################################################
############     ###  ######   ##        #####  #######  ##          ######################################################
############  ##  ##  #####  #  ####  #########  #####  #######  ##########################################################
############     ###  ####  ###  ###  ##########  ###  ###  ###  ##########################################################
############  ######  #####  #  ####  ###########  #  #####  #  ###########################################################
############  ######     ###   #####  ############   #######   ############################################################
###########################################################################################################################

    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
        
    if chain=="alpha":
        
	if tag_set == "original":       
	  tags_v = open("tags_trav.txt", "rU")
	  tags_j = open("tags_traj.txt", "rU")
 	elif tag_set == "extended":       
	  tags_v = open("exttags_trav.txt", "rU")
	  tags_j = open("exttags_traj.txt", "rU")       

        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v,num_j))
        for line in handle:
	  
            elements = line.rstrip("\n")
            
            if (tag_set == "original" and int(elements.split(',')[0]) >= 46) or (tag_set == "original" and int(elements.split(',')[1]) >= 49):
	      continue            

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))

	if tag_set == "original":
	  gene_list_v = ('V1-1','V1-2','V10','V12-1','V12-2','V12-3','V13-1','V13-2','V14/D4','V16','V17','V18','V19','V2','V20','V21','V22','V23/D6','V24','V25','V26-1','V26-2','V27','V29/DV5','V3','V30','V34','V35','V36/DV7','V38-1','V38-2/DV8','V39','V4','V40','V41','V5','V6','V7','V8-1','V8-2/8-4','V8-3','V8-6','V9-1','V9-2','DV1','DV2','DV3')
	  gene_list_j = ('J10','J11','J12','J13','J14','J15','J16','J17','J18','J20','J21','J22','J23','J24','J26','J27','J28','J29','J3','J30','J31','J32','J33','J34','J36','J37','J38','J39','J4','J40','J41','J42','J43','J44','J45','J46','J47','J48','J49','J5','J50','J52','J53','J54','J56','J57','J6','J7','J8','J9')
	  
	elif tag_set == "extended":
	  
	  gene_list_v = ('V1-1', 'V1-2', 'V10', 'V12-1', 'V12-2', 'V12-3', 'V13-1', 'V13-2', 'V14/DV4', 'V16', 'V17', 'V18', 'V19', 'V2', 'V20', 'V21', 'V22', 'V23/DV6', 'V24', 'V25', 'V26-1', 'V26-2', 'V27', 'V29/DV5', 'V3', 'V30', 'V34', 'V35', 'V36/DV7', 'V38-1', 'V38-2/DV8', 'V39', 'V4', 'V40', 'V41', 'V5', 'V6', 'V7', 'V8-1', 'V8-2/8-4', 'V8-3', 'V8-6', 'V9-1', 'V9-2', 'DV1', 'DV2', 'DV3', 'V8-7', 'V8-5', 'V11', 'V15', 'V28', 'V31', 'V32', 'V33', 'V37')
	  chromosome_order = [0, 1, 13, 24, 32, 35, 36, 37, 38, 42, 2, 49, 3, 39, 40, 6, 4, 48, 7, 8, 43, 50, 5, 41, 9, 10, 11, 12, 14, 15, 16, 17, 44, 18, 19, 20, 47, 22, 51, 23, 25, 52, 53, 54, 21, 26, 27, 28, 55, 29, 30, 31, 33, 34, 45, 46]
          gene_list_v = [ gene_list_v[i] for i in chromosome_order]
          
          gene_list_j = ('J10', 'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J20', 'J21', 'J22', 'J23', 'J24', 'J26', 'J27', 'J28', 'J29', 'J3', 'J30', 'J31', 'J32', 'J33', 'J34', 'J36', 'J37', 'J38', 'J39', 'J4', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J5', 'J50', 'J52', 'J53', 'J54', 'J56', 'J57', 'J6', 'J7', 'J8', 'J9', 'J1', 'J2', 'J19', 'J25', 'J35', 'J51', 'J55', 'J58', 'J59', 'J60', 'J61')
	  chromosome_order = [60, 59, 58, 57, 45, 44, 56, 43, 42, 41, 55, 40, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 27, 26, 25, 24, 54, 23, 22, 21, 20, 19, 17, 16, 15, 14, 53, 13, 12, 11, 10, 9, 52, 8, 7, 6, 5, 4, 3, 2, 1, 0, 49, 48, 47, 46, 39, 28, 18, 51, 50]
	  gene_list_j = [ gene_list_j[i] for i in chromosome_order]
          
            
        pos_v = np.arange(num_v)+ 1
        pos_j = np.arange(num_j)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_v)
        plt.xticks( pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, rotation='vertical', fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":

	if tag_set == "original":
	  tags_v = open("tags_trbv.txt", "rU")
	  tags_j = open("tags_trbj.txt", "rU")
 	elif tag_set == "extended":
	  tags_v = open("exttags_trbv.txt", "rU")
	  tags_j = open("exttags_trbj.txt", "rU")
	  
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v,num_j))
        for line in handle:

            elements = line.rstrip("\n")

            if (tag_set == "original" and int(elements.split(',')[0]) >= 45) or (tag_set == "original" and int(elements.split(',')[1]) >= 12):
	      continue	  

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))

	if tag_set == "original":        
	  gene_list_v = ('V10-1','V10-2','V10-3','V11-1','V11-2','V11-3','V12-3/V12-4','V12-5','V13','V14','V15','V16','V18','V19','V2','V20-1','V24-1','V25-1','V27','V28','V29-1','V3-1','V30','V4-1','V4-2','V4-3','V5-1','V5-4','V5-5','V5-6','V5-8','V6-1','V6-4','V6-5','V6-6','V6-8','V6-9','V7-2','V7-3','V7-4','V7-6','V7-7','V7-8','V7-9','V9')
	  gene_list_j = ('J1-1','J1-2','J1-3','J1-4','J1-5','J1-6','J2-1','J2-2','J2-3','J2-4','J2-5','J2-6','J2-7')
	  
	elif tag_set == "extended":

	  gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27', 'V28', 'V29-1', 'V3-1', 'V30', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9', 'V17', 'V23-1', 'V5-3', 'V5-7', 'V6-7', 'V7-1', 'V1', 'V12-1', 'V12-2', 'V21-1', 'V22-1', 'V26', 'V5-2', 'V7-5', 'V8-1', 'V8-2', 'VA', 'VB')
	  chromosome_order = [0, 1, 13, 24, 32, 35, 36, 37, 38, 42, 2, 49, 3, 39, 40, 6, 4, 48, 7, 8, 43, 50, 5, 41, 9, 10, 11, 12, 14, 15, 16, 17, 44, 18, 19, 20, 47, 22, 51, 23, 25, 52, 53, 54, 21, 26, 27, 28, 55, 29, 30, 31, 33, 34, 45, 46]
          gene_list_v = [ gene_list_v[i] for i in chromosome_order]
          
          gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7', 'J2-2P')
	  chromosome_order = [0, 1, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11, 12]
	  gene_list_j = [ gene_list_j[i] for i in chromosome_order]
          
          
          
          
        pos_v = np.arange(num_v)+ 1
        pos_j = np.arange(num_j)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_v)
        plt.xticks( pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags_v.close()
    tags_j.close()

def plot_insert_lengths( handle, savefilename="InsertLengths" ):

###########################################################################################################################
############     ###  ######   ##        ####    ##    ####  ##    ########################################################
############  ##  ##  #####  #  ####  ########  ###  #  ###  ###  #########################################################
############     ###  ####  ###  ###  ########  ###  ##  ##  ###  #########################################################
############  ######  #####  #  ####  ########  ###  ###  #  ##############################################################
############  ######     ###   #####  #######    ##  ####    ##############################################################
###########################################################################################################################
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    maxi = 500
    insert_lengths = [0]*maxi

    for line in handle:
        elements = line.rstrip("\n")

        classifier = elements.split(',')
        if len(classifier) == 6:
            insert_lengths[len(classifier[4].replace(' ',''))] += 1
        else:
            insert_lengths[0] += 1

    total = sum(insert_lengths)
    for i in range(len(insert_lengths)):
        insert_lengths[i] = insert_lengths[i] / float(total)

    ind = np.arange(maxi)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, insert_lengths, width, color='blue')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of nucleotides', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.xlim((0,50))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    

##################
#### ANALYSIS ####
##################

if (len(sys.argv) == 3):

  # chain explicitly provided
  
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
  
  # chain read from input file name
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
  print "e.g.: python PlotDCR.py FILE.txt a"
  print "or: python DCR.py FILE_beta.freq"       

name_results = filename.split(".")[0] 

currentpath = os.getcwd()

tag_set = default_tags

if platform.system() == 'Windows':
    newpath = currentpath+'\\Plots_'+str(name_results)+'\\' ## Ensure correct for specified platform
    if not os.path.exists(newpath):
        os.makedirs(newpath)
elif platform.system() == 'Linux':
    newpath = currentpath+'/Plots_'+str(name_results)+'/' ## Ensure correct for specified platform
    if not os.path.exists(newpath):
        os.makedirs(newpath)
elif platform.system() == 'Darwin':
    newpath = currentpath+'/Plots_'+str(name_results)+'/' ## Ensure correct for specified platform
    if not os.path.exists(newpath):
        os.makedirs(newpath)

print 'Plotting results of the analysis...'

plot_v_usage(open(filename, "rU"), chain=str(chain), savefilename = newpath+str(name_results)+'_Vusage', order="chromosome")
plot_j_usage(open(filename, "rU"), chain=str(chain), savefilename = newpath+str(name_results)+'_Jusage', order="chromosome")
plot_del_v(open(filename, "rU"), savefilename = newpath+str(name_results)+'_Vdels')
plot_del_j(open(filename, "rU"), savefilename = newpath+str(name_results)+'_Jdels')
plot_vj_joint_dist(open(filename, "rU"), chain=str(chain), savefilename = newpath+str(name_results)+'_VJusage')
plot_insert_lengths(open(filename, "rU"), savefilename = newpath+str(name_results)+'_InsertLengths')

print 'All plots (using', tag_set, 'tags) successfully saved to', newpath

