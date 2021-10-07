# This program will return a list of Tajima's D and a list of position for a given input VCF file with a window size
import argparse
import numpy as np
import pickle
import vcf
import warnings
warnings.filterwarnings('ignore')


# Get sample size, # of segregating site, base pair positions, and nucleotide diversities
def get_info(file_name):
    base_pair_positions = []
    nucleotide_diversities = []
    reader = vcf.Reader(open(file_name, 'r'))
    sample_size = len(reader.samples)*2
    for record in reader:
        base_pair_positions.append(record.POS)
        nucleotide_diversities.append(record.nucl_diversity)
    segregating_sites = len(base_pair_positions)
    return [sample_size, segregating_sites, base_pair_positions, nucleotide_diversities]


# Prepare constants for calculating the Tajima's D (n is the sample size)
def prepare_tajimas_d(n):
    a_1 = 0
    a_2 = 0
    for i in range(1, n):
        a_1 += 1/i
        a_2 += 1/i**2
    b_1 = (n + 1)/(3*(n - 1))
    b_2 = 2*(n**2 + n + 3)/(9*n*(n - 1))
    c_1 = b_1 - 1/a_1
    c_2 = b_2 - (n + 2)/(a_1*n) + a_2/a_1**2
    e_1 = c_1/a_1
    e_2 = c_2/(a_1**2 + a_2)
    return [a_1, e_1, e_2]


# Calculate the Tajima's D ( k is the nucleotide diversity and s is the # of segregating site)
def get_tajimas_d(k, s, a_1, e_1, e_2):
    if (np.sqrt(e_1*s + e_2*s*(s - 1))) == 0:
        return float('nan')
    else:
        return (k - s/a_1)/np.sqrt(e_1*s + e_2*s*(s - 1))


# Calculate the Tajima's Ds for parsed sequences as well as corresponding base pair positions
def analyze_parsed_sequence(sample_size, segregating_sites, base_pair_positions, nucleotide_diversities, window_size):
    [a_1, e_1, e_2] = prepare_tajimas_d(sample_size)
    parsed_positions = []
    tajimas_ds = []
    num = segregating_sites//window_size
    for index in range(num):
        parsed_positions.append(np.average(base_pair_positions[index*window_size:(index + 1)*window_size]))
        k = np.sum(nucleotide_diversities[index*window_size:(index + 1)*window_size])
        tajimas_ds.append(get_tajimas_d(k, window_size, a_1, e_1, e_2))
    if segregating_sites % window_size != 0:
        parsed_positions.append(np.average(base_pair_positions[num*window_size:]))
        k = np.sum(nucleotide_diversities[num*window_size:])
        tajimas_ds.append(get_tajimas_d(k, segregating_sites % window_size, a_1, e_1, e_2))
    return [parsed_positions, tajimas_ds]


# The main function.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Balancing selection analysis.')
    parser.add_argument(dest='file_name', help='Please enter the VCF filename.')
    parser.add_argument(dest='window_size', help='Please enter the window size.', type=int)
    args = parser.parse_args()
    print('The input file is: ' + str(args.file_name))
    print('The input window size is: ' + str(args.window_size))

    [Sample_size, Segregating_sites, Base_pair_positions, Nucleotide_diversities] = get_info(args.file_name)
    print('The sample size is: ' + str(Sample_size))
    print('The number of segregating site is: ' + str(Segregating_sites))
    print('The true length of the scaffold is: ' + str(Base_pair_positions[-1]))

    [Parsed_positions, Tajimas_ds] = analyze_parsed_sequence(
        Sample_size, Segregating_sites, Base_pair_positions, Nucleotide_diversities, args.window_size)
    with open(f'{args.file_name}_{args.window_size}_scores.pkl', 'wb') as t:
        pickle.dump(Tajimas_ds, t)
    with open(f'{args.file_name}_{args.window_size}_positions.pkl', 'wb') as p:
        pickle.dump(Parsed_positions, p)
    print('The number of parsed fragment is: ' + str(len(Tajimas_ds)))
