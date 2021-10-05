# Updated Tajima's D analysis
import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pprint as ppt
import vcf
import pickle
import os
warnings.filterwarnings("ignore")


# Get sample size, base pair positions, and allele frequencies from the VCF file
def get_info(file_name):
    base_pair_positions = []
    allele_frequencies = []
    reader = vcf.Reader(open(file_name, 'r'))
    sample_size = len(reader.samples)
    for record in reader:
        base_pair_positions.append(record.POS)
        if len(record.INFO["AF"]) != 1:
            print("Non-binary polymorphism detected at position: " + str(record.POS))
            allele_frequencies.append(1.0)
        else:
            allele_frequencies.append(record.INFO["AF"][0])
    return [sample_size, base_pair_positions, allele_frequencies]


# Calculate nCr instead of using the math package
def comb(n, k):
    ntok = 1
    ktok = 1
    for t in range(1, min(k, n - k) + 1):
        ntok *= n
        ktok *= t
        n -= 1
    return ntok//ktok


# Calculate the nucleotide diversity using allele frequencies
def get_nucleotide_diversity(allele_frequencies, sample_size):
    total_difference = 0
    for frequency in allele_frequencies:
        total_difference += frequency*(1 - frequency)*(sample_size**2)
    return total_difference/comb(sample_size, 2)


# Prepare constants for calculating the Tajima's D where n is the sample size
def prepare_tajimas_d(n):
    a_1 = 0
    a_2 = 0
    for number in range(1, n):
        a_1 += 1/number
        a_2 += 1/number**2
    b_1 = (n + 1)/(3*(n - 1))
    b_2 = 2*(n**2 + n + 3)/(9*n*(n - 1))
    c_1 = b_1 - 1/a_1
    c_2 = b_2 - (n + 2)/(a_1*n) + a_2/a_1**2
    e_1 = c_1/a_1
    e_2 = c_2/(a_1**2 + a_2)
    return [a_1, e_1, e_2]


# Calculate the Tajima's D where k is the nucleotide diversity and s is the number of segregating site
def get_tajimas_d(k, s, a_1, e_1, e_2):
    if (np.sqrt(e_1*s + e_2*s*(s - 1))) == 0:
        return float("nan")
    else:
        return (k - s/a_1)/np.sqrt(e_1*s + e_2*s*(s - 1))


def parse_the_frequency(allele_frequencies, window_size):
    num = len(allele_frequencies)//window_size
    parsed_frequencies = [None]*num
    for index in range(num):
        parsed_frequencies[index] = allele_frequencies[index*window_size:(index + 1)*window_size]
    if len(allele_frequencies) % window_size != 0:
        parsed_frequencies.append(allele_frequencies[(index + 1)*window_size:])
    return parsed_frequencies


# The main function.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Balancing selection detector.')
    parser.add_argument(dest='file_name', help='Please enter the filename.')
    parser.add_argument(dest='window_size', help='Please enter the window size.', type=int)
    args = parser.parse_args()

    print("The input file is: " + str(args.file_name))

    [Sample_size, Base_pair_positions, Allele_frequencies] = get_info(args.file_name)
    print("The sample size is: " + str(Sample_size))
    print("The number of polymorphic base pair is: " + str(len(Base_pair_positions)))
    print("The true length of the genome is: " + str(Base_pair_positions[-1]))

    Pi_value = get_nucleotide_diversity(Allele_frequencies, Sample_size)
    print("The nucleotide diversity is: " + str(Pi_value))

    [A_1, E_1, E_2] = prepare_tajimas_d(Sample_size)
    Tajima_D = get_tajimas_d(Pi_value, len(Base_pair_positions), A_1, E_1, E_2)
    print("The Tajima's D of the genome is: " + str(Tajima_D))

    Parsed_frequencies = parse_the_frequency(Allele_frequencies, args.window_size)
    print("The number of row is: " + str(len(Parsed_frequencies)))
    print("The number of column is: " + str(len(Parsed_frequencies[0])))
    ppt.pprint(Parsed_frequencies[-1])
