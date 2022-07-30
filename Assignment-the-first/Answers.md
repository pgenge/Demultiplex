# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 | 101 | Phred 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | Phred 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | Phred 33 |
| 1294_S1_L008_R4_001.fastq.gz | Read 3 | 101 | Phred 33 |


2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. ![Index 1](https://github.com/pgenge/Demultiplex/blob/master/Assignment-the-first/index1_quality.png)
    3. ![Index 2](https://github.com/pgenge/Demultiplex/blob/master/Assignment-the-first/index2_quality.png)
    4. ![Read 1](https://github.com/pgenge/Demultiplex/blob/master/Assignment-the-first/read1_quality.png)
    5. ![Read 2](https://github.com/pgenge/Demultiplex/blob/master/Assignment-the-first/read2_quality.png)
    
## Part 2
1. Define the problem
`We have 4 files: R1, R2, R3, R4, and we need to separate the records in the files based on index (i.e.each sample). We also need to separate our reads based on quality criteria. Our problem is that some index sequences in R2/R3 don't match our known index sequences, some contain low quality base calls, some contain unknown base calls (N), and some index1-index2 index matches for each record do not exist (index hopping has occurred). We want to report high quality data and to do so we have to pull out all of these things and only report what meets our criteria of high quality index reads (Q>/=30), matched to our known indexes/don't contain unknown basecalls, and index pairs match each other for each indexed sample.`

2. Describe output
```
Total number of files = 52 (or more if we want to pull out samples that have medium quality reads or 1 basecall unknown but minimum is 52)
Each insert sequence read from R1/R4 correspoding to the index reads from R2/R3 will output the following:
    R1 matched
    R2 matched
        Total # of files = 48 (2 for each sample/known index)
Each matched file will contain reads for that index/sample that we looked for based on our known index list that have index1-index2 match and high quality reads for that index with no unknown basecalls.
If the index reads in our R2/R3 files do not have corresponding index1-index2 matches we'll write out the R1/R4 insert seq read/record respectively to:
    R1 Unmatched
    R2 Unmatched
Furthermore if the index reads in R2/R3 do not match our known indexes, contain unknown basecalls, or low quality basecalls, we'll write out the R1/R4 insert seq read/record respectively to:
    R1 Unknown
    R2 Unknown

```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
```
Make 1 set with forward seqeuence of 24 known indexes 
write a function to do reverse complement

Open all 4 fastq files & read in parallel
for nth record (meaning every record in file): 
        check if index 1 in set
            if NOT in set
                write R1/R4 record to unknown files R1_unknown R2_unknown
            if index 1 IS IN set
                check if revcomp(index 2) is also in set
                    if NOT send to unknown files corresponding to read
                    if index 2 is in set
                        check line 4 of R2/R3 records for quality score of index for Q>/=30
                            if below cutoff write to unknown files corresponding to read
                            if Q>30 then:
                                check if index 1 matches revcomp(index 2)
                                    if it's a match
                                        send R1/R4 record to matched files R1_matched R2_matched
                                    if NOT matched
                                        write R1/R4 record to unmatched files R1_unmatched R2_unmatched
                                            each time you a write out a record append index1_revcomp(index2) to each header
```
5. High level functions. For each function, be sure to include:
    1. Description/doc string

    2. Function headers (name and parameters)

    3. Test examples for individual functions

    4. Return statement
```
def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

Test example:
E = 33
A = 28
# = 2
convert_phred(E) == 33 etc.. 

def qual_score(phred_score: str) -> float:
    """This function calculates the average quality score of the whole phred string"""
    sum = 0
    for x in phred_score:
        sum += convert_phred(x)
    return sum/len(phred_score)


EA# 
```
