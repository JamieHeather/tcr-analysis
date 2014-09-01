#! /bin/bash

# Uniqify.sh

### PURPOSE ###

# Take a non-collapsed, decombined file and collapse it down into its unique assignations, with frequency scores
# Allows comparison of frequency metrics between pre- and post-processed files, producing files in same format that say CollapseTCRs.py produces

### INPUT/RUNNING ###

# Raw output from any version of decombinator, e.g. *.n12 files produced by vDCR

# Run:
# sh Uniqify.sh FILENAME.fq
# for x in *n12; do echo $x; sh Uniqify.sh $x; done

### OUTPUT ### 

# File containing 5 part identifier, with a sixth field of frequency of occurrence
# Output filename = FILENAME.uniq

in=$1
out=${1%.*}

awk 'BEGIN { FS = ", "} ; { print $1,$2,$3,$4,$5}' $in | sort | uniq -c | sort -rn | sed 's/^[ \t]*//' | awk '{print $2,$3,$4,$5,$6,$1}' | sed 's/^[ \t]*//;s/ /, /g' > $out.uniq




