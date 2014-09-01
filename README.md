#JamieHeather/tcr-analysis README.md v1

University College London, UK, 2014

This repository has been built to contain the suite of inter-related T-cell receptor (TCR) analysis scripts that I've developed in the course of my PhD, here made available as an adjunct to my doctoral thesis.

A number of these scripts were also used in the analysis of the data that made up the results of my paper (Heather et al, ***ADD REFERENCE***).

***ADD LINKS TO SEQUENCE READS/DATA***

Please note that a number of these scripts are either modifications to existing code (such as vDCR.py, which is my re-working of code produced by Dr. Niclas Thomas) or in close colloboration with other lab members, particularly Katharine Best.

===== STANDARD PIPELINE =====

##0) Deep-sequence TCR transcripts!

- See http://dx.doi.org/10.6084/m9.figshare.950994 for an overview of our amplification procedure

##1) PREPARE FASTQs 
###	• AddAllR2Hex.py
- Used to generate a barcoded file, where the first twelve bases represent a unique barcode relating to original cDNA molecule
- Transfers the first 6 bases of a given R2 fastq to the beginning of a R1 file, producing a third fastq
- For use if you only have the two demultiplexed paired end fastqs from the MiSeq
- Outputs a third '.fq' file

###    • DualIndexDemultiplexing.py
- Performs a similar role as AddAllR2Hex.py, but used when you additionally have the index read fastq
- Allows manual demultiplexing of samples along with simultaneous barcode preparation
- For use in scenarios where you don't trust Illumina's inbuilt demultiplexing
- Either requires manual reconstitution of fastqs from the bcl files, or to set the necessary MiSeq flags before running

 
##2) DECOMBINE 
###    • vDCR.py
- Built on standard Decombinator (DCR) v1.4
  - See https://github.com/uclinfectionimmunity/Decombinator/ and http://dx.doi.org/10.1093/bioinformatics/btt004
- Outputs additional fields beyond the typical 5:
  - ID of read
  - TCR nucleotide sequence (start of V tag to end of J tag)
  - TCR quality (same sequence window)
  - Barcode sequence
  - Barcode quality
- Outputs a '.n12' output file (to distinguish from output of regular dcr)
  - This is the required input format for error- and frequency-correction by CollapseTCRs.py
    
    
##3) ERROR-CORRECT
###    • CollapseTCRs.py
- Takes n12 files and outputs '.freq' files of unique TCRs, with frequencies
- First collapses TCRs by nucleotide sequence, correcting PCR and sequencing errors
- Next clusters barcodes to mitigate for PCR amplification, providing a more accurate frequency value
- NB various filters and thresholds can be tweaked to alter stringency of output results
    

##4) TRANSLATE
###    • CDR3ulator.py
- Generates CDR3 sequences of productive rearrangements
- Takes any Decombined output (.txt/.n12/.freq etc) as input
- Recapitulates entire nucleotide sequence, translates this and then extracts the CDR3 sequence
- Has two major options for controlling format of output
  - Users can include or ignore the final GXG motif 
  - Output can be just unique CDR3s (with their frequencies), or can be unique DCR assignations with their CDR3s/frequencies
    - Using this option means that CDR3s can appear on multiple lines, as they can be encoded by multiple DCRs
    
NOTES ON USAGE
- All scripts assume that any Decombined file has the TCR 5-part identifier as the first 5 fields, in the same order
  - Hence any downstream analysis must not alter these, only append to them
- Filenames are often used to provide information to downstream scripts, therefore:
  - try to maintain extensions where possible
  - don't remove "alpha" or "beta" from filenames 
  - don't introduce extra full stop characters (e.g. 'file.first-try.txt' would cause problems)
- All python scripts were written in and tested on Python version 2.7.6
- All R scripts were written in and tested on R version 3.0.2
 
 
====== OPTIONAL EXTRAS ======

##PLOTTING
###    • PlotDCR.py 
- Provides a means to obtain the V/J/pairing/deletion frequency plots produced by Decombinator v1.4 from any Decombinator output
    

##FORMATTING
###    • Uniqify.sh 
- Collapses any redundant decombined file (i.e. raw output of vDCR or Decombinator) into unique lines
- Collapses based purely on DCR 5-part identifier, with no error-correction 
- Can be used to reduce file size and find the number of unique TCR assignations in raw data     




