# Author: Palak Genge pgenge@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.8"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "ACTG"
RNA_bases = "ACUG"

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter)-33

def qual_score(phred_score):
    """This function calculates the average quality score of the whole phred string"""
    sum = 0
    for x in phred_score:
        sum += convert_phred(x)
    return sum/len(phred_score)

# def qual_score(phred_string: str) -> float:
#     """This function calculates the average quality score of the whole phred string"""
#     sum = 0
#     for x in phred_string:
#         sum += convert_phred(x)
#     avg_phred_score = sum/len(phred_string)
#     return avg_phred_score


DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)


# assert True == False

# assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
# assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
# assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
# assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
# print("Passed DNA and RNA tests")

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)




def oneline_fasta(line):
    '''This function concatenates the sequence lines of each record to be on the same line'''
    i = 0 #don't get the 1st line
    for line in fa2: #loop through lines in fa
        if '>' in line and i!=0: #if it is header and not the first line, write a new line so that you don't have just one line but header then seq
            return('\n') #write out the new line after header
        if '>' not in line: #if not header strip the new line characters
            line = line.strip('\n')
            return(line) #write out the seq line
        i += 1



#unit tests
if "bioinfo.py"  == "__main__":
    pass
    assert gc_content("GCGCGC") == 1, "messed up calc when all GC"
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
print("correctly calculated GC content")