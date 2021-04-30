# This file contains the calculation of nucleotide diversity.
import math

import numpy as np
import warnings
import argparse
from Parse_the_input import *
from pprint import pprint
warnings.filterwarnings("ignore")


def parse_to_allele(sequence, window_size):  # Parse each sequence into pieces with fixed length
    allele = []
    for i, sequence in enumerate(sequence):
        num = len(sequence)//window_size
        allele.append([]*num)
        for index in range(num):
            allele[i].append(sequence[index*window_size:(index + 1)*window_size])
        allele[i].append(sequence[num*window_size:])
    return allele


def pairwise_difference(first_sequence, second_sequence, threshold):  # Single nucleotide polymorphism
    count = np.sum(nucleotide_A != nucleotide_B for nucleotide_A, nucleotide_B
                   in zip(first_sequence, second_sequence))\
            + abs(len(first_sequence) - len(second_sequence))
    if count < threshold:
        count = 0
    return count


def nucleotide_diversity(sequences, threshold):  # The nucleotide diversity of a population
    pi_value = 0
    n = len(sequences)
    for first_sequence in sequences:
        for second_sequence in sequences:
            pi_value += pairwise_difference(first_sequence, second_sequence, threshold)
    pi_value = pi_value/math.comb(n, 2)
    return pi_value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genome reader.')
    parser.add_argument(dest='file_name', help='Please input the genome filename.')
    args = parser.parse_args()
    Sequences = file_parser(args.file_name)
    pprint(Sequences)
    Allele = parse_to_allele(Sequences, 50)
    pprint(Allele)
    Pi_value = nucleotide_diversity(Sequences, 100)
    pprint(Pi_value)
