# This file contains the calculation of nucleotide diversity

# Import Packages
from pprint import pprint
from numpy import *

# Import the FASTA genome file
file_name = "test1.fas"

# Convert the FASTA genome file into sequences
def read_file(genome):
    sequence = []
    with open(genome) as f:
        current_sequence = ""
        for line in f:
            if line[0] != '>':
                current_sequence = current_sequence + line.strip()
            else:
                if current_sequence:
                    sequence.append(current_sequence)
                    current_sequence = ""
    return sequence

# Parse full sequences into alleles with corresponding window size
def parse_to_allele(sequence,window_size):
    allele = []
    for i,sequence in enumerate(sequence):
        num = len(sequence)//window_size
        allele.append([]*num)
        for index in range(num):
            allele[i].append(sequence[index*window_size:(index + 1)*window_size])
        allele[i].append(sequence[num*window_size:])
    return allele

if __name__ == "__main__":
    Sequence = read_file(file_name)
    Allele = parse_to_allele(Sequence,60)
    pprint(Allele[2][3])