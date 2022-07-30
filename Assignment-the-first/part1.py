#!/usr/bin/env python

#import libs
import bioinfo
import numpy as np
import gzip
import argparse
import matplotlib.pyplot as plt


#argparse for calling options

def get_args():
    parser = argparse.ArgumentParser(description="get distribution of phred scores")
    parser.add_argument('-f', '--file', help = 'File name', type = str, required = True)
    parser.add_argument('-o', '--output', help = 'Output File Name For Histogram', required = True)
    parser.add_argument('-l', '--readlength', help = 'Read Length of the Sequence', type = int)
    parser.add_argument('-t', '--title', help = 'Title for Histogram')
    return parser.parse_args()

#allows method calling
args = get_args() 

#initialize list
qual_scores = [0] * args.readlength

#make a list with the running sum
def populate_list(file):
    """This function populates an empty list with the converted Phred score 
    and adds them together while indexing the position of the nucleotides"""
    linecount = 0 #start indexing at 0
    with gzip.open(args.file, "rt") as fastq: #open the file
        for line in fastq: #count lines
            if linecount%4 == 3: #checking if the remainder is 3
                line = line.strip('\n') 
                for j, x in enumerate(line): #count positions of nts
                    qual_scores[j] += bioinfo.convert_phred(x)  #
            linecount += 1    #increment 1 each time    
    return (qual_scores, linecount)

#call populate_list and pass it a file
qual_scores, num_lines = populate_list(args.file)
num_records = num_lines/4

#get average
for i, x in enumerate(qual_scores): #keep track of count, each value in the count pass my_list
   qual_scores[i] = x/num_records #get average for each position

#plot and save as png

x = range(0, args.readlength)
y = qual_scores

plt.subplots(1, figsize=(15, 12))
fig, ax = plt.subplots()
ax.bar(x, y)
ax.set_xlabel("Base Position in Read")
ax.set_ylabel("Mean Quality Score")
ax.set_title(f"Per Base Mean Quality Score Distribution\n{args.title}")

plt.savefig(args.output) #.png





