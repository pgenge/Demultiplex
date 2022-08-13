#!/usr/bin/env python

#file path to test files: 
# /projects/bgmp/pgenge/bioinfo/Bi622/Demultiplex/Assignment-the-first/R1.fq
# /projects/bgmp/pgenge/bioinfo/Bi622/Demultiplex/Assignment-the-first/R2.fq
# /projects/bgmp/pgenge/bioinfo/Bi622/Demultiplex/Assignment-the-first/R3.fq
# /projects/bgmp/pgenge/bioinfo/Bi622/Demultiplex/Assignment-the-first/R4.fq

#test command: ./demux.py R1.fq.gz R2.fq.gz R3.fq.gz R4.fq.gz indextest.txt

#import all required modules
import bioinfo
import numpy as np #could use but not necessary--I'm not doing arrays
import gzip
import argparse
import matplotlib.pyplot as plt

def get_args():
    #positional arguments--see demux_nb.txt
    parser = argparse.ArgumentParser(description="demux")
    parser.add_argument('r1_filename', default=None, type=str, help='Specify R1 file from sequencing run file.')
    parser.add_argument('r2_filename', default=None, type=str, help='Specify R2 file from sequencing run file.')
    parser.add_argument('r3_filename', default=None, type=str, help='Specify R3 file from sequencing run file.')
    parser.add_argument('r4_filename', default=None, type=str, help='Specify R4 file from sequencing run file.')
    parser.add_argument('indexfile', default=None, type=str, help='Specify known index file.')
    #optional argument
    parser.add_argument('-o', '--outputdir', help = 'Specify output directory for files')
    return parser.parse_args()

#allows method calling
args = get_args()

#open all input files in parallel all at once
indexfile = open(args.indexfile, "r")
r1 = gzip.open(args.r1_filename, "rt")
r2 = gzip.open(args.r2_filename, "rt")
r3 = gzip.open(args.r3_filename, "rt")
r4 = gzip.open(args.r4_filename, "rt")

#write reverse complement function
def revcomp(DNA: str):
    """This function takes a DNA sequence and returns the reverse complement"""
    revcomp_seq = ""
    bases_dict = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"} 
    for base in DNA:
        revcomp_seq += bases_dict[base]
    return (revcomp_seq[::-1]) #this allows to read the string backwards

#make set with indexes.txt
indexes = set()
#make reverse comp set with indexes.txt
index_revcomp = set()
indexfile.readline()
for line in indexfile:
    #get rid of new line
    line = line.strip("\n")
    #split the line by tabs
    line = line.split("\t")
    #grab the 5th line (index seq)
    indexes.add(line[4])
    index_revcomp.add(revcomp(line[4]))

#define function for getting each record
def get_records(filehandle):
    """This function gets each line of a fastq record and sets it equal to a variable, 
    categorizes them based on content and returns each line as separate variable"""
    #set each line to a variable
    header = filehandle.readline().strip("\n")
    seq = filehandle.readline().strip("\n")
    plus = filehandle.readline().strip("\n")
    qscore = filehandle.readline().strip("\n")
    return header, seq, plus, qscore 

#make diction with types of non-matched demux files
category = {"unknown","hopped"}

#writing out the files with appropriate file handles
#make an empty dictionary to populate below
write_files = {}
#this is for writing out to the unknown and hopped files
for type in category:
    r1_filename = open(f"{args.outputdir}/{type}_R1.fq", "w")
    r2_filename = open(f"{args.outputdir}/{type}_R2.fq", "w")
    write_files[type] = [r1_filename, r2_filename]
#make dictionary with known seq of index as key, and the value is the filehandle
#need a dict with the key index and the value as the open statement to open and write files to index match file
for index in indexes:
    r1_filename = open(f"{args.outputdir}/{index}_{index}_R1.fq", "w")
    r2_filename = open(f"{args.outputdir}/{index}_{index}_R2.fq", "w")
    write_files[index] = [r1_filename, r2_filename]


#get stats initialize these at zero to capture totals while looping
totalnumrecords = 0
unknown = 0
unmatched = 0
matched = 0

# make empty dict for number of matched reads and unmatched reads
nummatchedreads = {}
numunmatchedreads = {}


#loop through each file and filter for indexes based on conditions
while True:
    #get records for every file
    read1_head, read1_seq, read1_plus, read1_qs = get_records(r1)
    index1_head, index1_seq, index1_plus, index1_qs = get_records(r2)
    index2_head, index2_seq, index2_plus, index2_qs = get_records(r3)
    read2_head, read2_seq, read2_plus, read2_qs = get_records(r4)
    #BREAK so that while loop knows to end
    if read1_head == "":
        break
#count the number of total records before filtering
    totalnumrecords += 1
#creating a variable with an f-string for writing out file headers easily
    records_r1 = (f"{read1_head}\n{read1_seq}\n{read1_plus}\n{read1_qs}\n")
    records_r2 = (f"{read2_head}\n{read2_seq}\n{read2_plus}\n{read2_qs}\n")
#checking if the index 1 (i7) in run file is a known index in given index file
    if index1_seq not in indexes:
        #if not in index file writing to unknown files for each R1-R2
        write_files["unknown"][0].write(records_r1)
        write_files["unknown"][1].write(records_r2)
        #counting unknowns
        unknown += 1
#checking if index 2 (i5) in run file is a known index in given index file    
    elif index2_seq not in index_revcomp:
        #writing to unknown file when it meets the above condition not in set of known index 2s
        write_files["unknown"][0].write(records_r1)
        write_files["unknown"][1].write(records_r2)
        #counting unknowns
        unknown += 1
#taking the average of the qual score of the index 1 read and checking if below 30 (cutoff)
    elif bioinfo.qual_score(index1_qs) < 30:
        #if it doesn't meet average > 30 criteria-write to unknown and increment counter
        write_files["unknown"][0].write(records_r1)
        write_files["unknown"][1].write(records_r2)
        unknown += 1
#taking the average of the qual score of the index 2 read and checking if below 30 (cutoff)
    elif bioinfo.qual_score(index2_qs) < 30:
        #if it doesn't meet average > 30 criteria-write to unknown and increment counter
        write_files["unknown"][0].write(records_r1)
        write_files["unknown"][1].write(records_r2)
        unknown += 1
#checking if index 1 does not match reverse complement index 2
    elif index1_seq != revcomp(index2_seq):
        #if above is true=the indexes do not match then writing to a "hopped" file
        write_files["hopped"][0].write(records_r1)
        write_files["hopped"][1].write(records_r2)
        #incrementing counter
        unmatched += 1
        #populating dictionary with index1_index2 seqs as the key and the counts as value
        if f"{index1_seq}_{revcomp(index2_seq)}" in numunmatchedreads:
            numunmatchedreads[f"{index1_seq}_{revcomp(index2_seq)}"] += 1
        else:
            numunmatchedreads[f"{index1_seq}_{revcomp(index2_seq)}"] = 1
#checking if the indexes are dual matched (index1 = reverse complement index 2)
    elif index1_seq == revcomp(index2_seq):
        #appending the "index1_revcomp(index2)" to the header for each record
        read1_head += f"_{index1_seq}_{revcomp(index2_seq)}"
        read2_head += f"_{index1_seq}_{revcomp(index2_seq)}"
        #new write files to include the appended header
        write_files[index1_seq][0].write(f"{read1_head}\n{read1_seq}\n{read1_plus}\n{read1_qs}\n")
        write_files[index1_seq][1].write(f"{read2_head}\n{read2_seq}\n{read2_plus}\n{read2_qs}\n")
        #incrementing counter for stats
        matched += 1
        #populating dictionary with just index 1 seq as key because index 1 = revcomp(index2) to count frequency
        if index1_seq in nummatchedreads:
            nummatchedreads[index1_seq] += 1
        else:
            nummatchedreads[index1_seq] = 1

#loop through files and close each one
for file in write_files:
    write_files[file][0].close()
    write_files[file][1].close()
    
#write the stats file and capture important information
with open(f"{args.outputdir}/run_stats.txt", "w") as statsfile:
    statsfile.write(f"Total Number of Reads:{totalnumrecords}\n")
    statsfile.write (f"Total Number of Dual Matched Reads: {matched} Percent of Dual-Matched Reads: {(matched/totalnumrecords)*100}%\n")
    statsfile.write(f"Total Number of Uknown Reads:{unknown} Percent of Uknown Reads: {(unknown/totalnumrecords)*100}%\n")
    statsfile.write(f"Total Number of Unmatched Read: {unmatched} Percent of Unmatched Reads: {(unmatched/totalnumrecords)*100}%\n")
    statsfile.write(f"Index_Seq\tNumber_of_Reads\tPercentage_of_Total_Reads\n")
    #writing out the dual matched indexes and frequency counts + percentage
    for index in nummatchedreads:
        statsfile.write(f"{index}\t{nummatchedreads[index]}\t{(nummatchedreads[index]/totalnumrecords)*100}%\n")
    statsfile.write(f"Index_Seq_Hopped\tNumber_of_Reads\tPercentage_of_Total_Reads\n")
    #writing out the UNmatched indexes and frequency counts + percentage
    for index in numunmatchedreads:
        statsfile.write(f"{index}\t{numunmatchedreads[index]}\t{(numunmatchedreads[index]/totalnumrecords)*100}%\n")

#plot percentage of reads as pie chart
labels = 'Dual Matched Reads', 'Unknown Reads','Unmatched Reads (Index Hopping)'
sizes = [(matched/totalnumrecords)*100, (unknown/totalnumrecords)*100,(unmatched/totalnumrecords)*100]

fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.savefig(f"{args.outputdir}/frequency_of_type_of_read.png")

#plot frequency of dual matched indexes
matchedindex_list = list(nummatchedreads.keys())
matchedcounts_list = list(nummatchedreads.values())
x = matchedindex_list
y = matchedcounts_list

fig, ax = plt.subplots()
ax.bar(x, y)
ax.set_xlabel("Indexes")
ax.set_ylabel("Counts (Number of Reads)")
ax.set_title("Number of Reads per dual-matched Index")
plt.xticks(rotation=50)
plt.tight_layout()

plt.savefig(f"{args.outputdir}/frequency_of_dual-matched_indexes.png") #.png
