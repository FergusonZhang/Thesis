# This file contains the calculation of nucleotide diversity
import numpy as np
import warnings
from pprint import pprint
from numpy import *
warnings.filterwarnings("ignore")

# Import the FASTA genome file
file_name = "Test1.fas"

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
        if current_sequence:
            sequence.append(current_sequence)
        return sequence

# Parse sequences into alleles with corresponding window size
def parse_to_allele(sequence,window_size):
    allele = []
    for i,sequence in enumerate(sequence):
        num = len(sequence)//window_size
        allele.append([]*num)
        for index in range(num):
            allele[i].append(sequence[index*window_size:(index + 1)*window_size])
        allele[i].append(sequence[num*window_size:])
    return allele

# Calculate the pairwise difference between two sequences with justification
def pairwise_difference(first_sequence,second_sequence,threshold):
    count = np.sum(x != y for x,y in zip(first_sequence,second_sequence))\
            + abs(len(first_sequence) - len(second_sequence))
    if count < threshold:
        count = 0
    return count

# Calculate the nucleotide diversity of a given allele set
def nucleotide_diversity(sequence,threshold):
    pi_value = 0
    # TODO: more implements of the algorithm
    frequncy = 1/len(sequence)
    for first_sequence in sequence:
        for second_sequence in sequence:
            pi_value += pairwise_difference(first_sequence,second_sequence,threshold)*frequncy
    return pi_value

if __name__ == "__main__":
    Sequence = read_file(file_name)
    pprint(Sequence)
    Allele = parse_to_allele(Sequence,50)
    pprint(Allele)
    Pi_value = nucleotide_diversity(Sequence,100)
    pprint(Pi_value)