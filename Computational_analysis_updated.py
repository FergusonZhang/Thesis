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


# Get sample size, segregating sites, base pair positions, and allele frequencies from the VCF file
def get_info(file_name):
    base_pair_positions = []
    allele_frequencies = []
    test_pi = 0
    reader = vcf.Reader(open(file_name, 'r'))
    sample_size = len(reader.samples)*2
    for record in reader:
        base_pair_positions.append(record.POS)
        test_pi += record.nucl_diversity
        if len(record.INFO["AF"]) != 1:
            print("Non-binary polymorphism detected at position: " + str(record.POS))
            allele_frequencies.append(1.0)
        else:
            allele_frequencies.append(record.INFO["AF"][0])
    segregating_sites = len(base_pair_positions)
    return [sample_size, segregating_sites, base_pair_positions, allele_frequencies, test_pi]


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


# Calculate the Tajima's D where k is the nucleotide diversity and s is the segregating sites
def get_tajimas_d(k, s, a_1, e_1, e_2):
    if (np.sqrt(e_1*s + e_2*s*(s - 1))) == 0:
        return float("nan")
    else:
        return (k - s/a_1)/np.sqrt(e_1*s + e_2*s*(s - 1))


# Parse the base pair positions and the frequencies using the input window size
# def parse_the_sequence(base_pair_positions, allele_frequencies, window_size):
#     num = len(allele_frequencies)//window_size
#     parsed_positions_raw = [None]*num
#     parsed_frequencies = [None]*num
#     for index in range(num):
#         parsed_frequencies[index] = allele_frequencies[index*window_size:(index + 1)*window_size]
#         parsed_positions_raw[index] = base_pair_positions[index*window_size:(index + 1)*window_size]
#     # if len(allele_frequencies) % window_size != 0:
#     #     parsed_frequencies.append(allele_frequencies[(index + 1)*window_size:])
#     #     parsed_positions.append(base_pair_positions[(index + 1)*window_size:])
#     parsed_positions = []
#     for position in parsed_positions_raw:
#         parsed_positions.append(np.average(position))
#     return [parsed_positions, parsed_frequencies]


# Calculate the Tajima's Ds for parsed sequence with a fixed window size
def analyze_parsed_sequence(sample_size, segregating_sites, base_pair_positions, allele_frequencies, window_size):
    [a_1, e_1, e_2] = prepare_tajimas_d(sample_size)
    num = segregating_sites//window_size
    tajima_ds = []
    parsed_positions = []
    for index in range(num):
        parsed_positions.append(np.average(base_pair_positions[index*window_size:(index + 1)*window_size]))
        k = get_nucleotide_diversity(allele_frequencies[index*window_size:(index + 1)*window_size], sample_size)
        tajima_ds.append(get_tajimas_d(k, window_size, a_1, e_1, e_2))
    return [parsed_positions, tajima_ds]


# The main function.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Balancing selection detector.')
    parser.add_argument(dest='file_name', help='Please enter the filename.')
    parser.add_argument(dest='window_size', help='Please enter the window size.', type=int)
    args = parser.parse_args()
    print("The input file is: " + str(args.file_name))
    print("The input window size is: " + str(args.window_size))

    [Sample_size, Segregating_sites, Base_pair_positions, Allele_frequencies, Test_pi] = get_info(args.file_name)
    print("The sample size is: " + str(Sample_size))
    print("The number of segregating site is: " + str(Segregating_sites))
    print("The true length of the genome is: " + str(Base_pair_positions[-1]))

    Pi_value = get_nucleotide_diversity(Allele_frequencies, Sample_size)
    print("The nucleotide diversity is: " + str(Pi_value))
    print("The alternative diversity is: " + str(Test_pi))
    #
    # [A_1, E_1, E_2] = prepare_tajimas_d(Sample_size)
    # Tajima_D = get_tajimas_d(Pi_value, len(Base_pair_positions), A_1, E_1, E_2)
    # print("The Tajima's D of the genome is: " + str(Tajima_D))

    # [Parsed_positions, Parsed_frequencies] = \
    #     parse_the_sequence(Base_pair_positions, Allele_frequencies, args.window_size)
    # print("The number of row is: " + str(len(Parsed_frequencies)))
    # print("The number of column is: " + str(len(Parsed_frequencies[0])))
    # print("The length of the last element is: " + str(len(Parsed_frequencies[-1])))

    [Parsed_positions, Tajima_Ds] = \
        analyze_parsed_sequence(Segregating_sites, Base_pair_positions, Allele_frequencies, args.window_size)
    print("The number of fragment is: " + str(len(Tajima_Ds)))

    plt.plot(Parsed_positions, Tajima_Ds, color='black', linestyle='dashed', linewidth=0.5,
             marker='.', markerfacecolor='blue', markersize=2.5)
    plt.xlabel("Base Pair Position")
    plt.ylabel("Tajima's D")
    plt.title(f"{args.file_name} Balancing Selection Analysis")
    plt.savefig(f"{args.file_name} Tajima's D.png")
