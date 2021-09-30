# Tajima's D analysis
import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pprint as ppt
import vcf
import pickle
import os
warnings.filterwarnings("ignore")


# Read sequences from different file types.
def read_sequences(file_name):
    sequences = []
    file_type = file_name.split(".")[1]
    if file_type == "fas":
        current_sequence = ""
        with open(file_name) as f:
            for line in f:
                if line[0] != ">":
                    current_sequence = current_sequence + line.strip()
                else:
                    if current_sequence:
                        sequences.append(current_sequence)
                        current_sequence = ""
            if current_sequence:
                sequences.append(current_sequence)
        return sequences
    elif file_type == "vcf":
        reader = vcf.Reader(open(file_name, 'r'))
        sequences = [None]*(len(reader.samples)*2)
        for record in reader:
            sample_index = 0
            for sample_name in reader.samples:
                allele = record.genotype(sample_name).gt_bases
                sequences[sample_index] = str(sequences[sample_index]) + allele[0]
                sequences[sample_index + 1] = str(sequences[sample_index + 1]) + allele[2]
                sample_index += 2
        return sequences
    else:
        print("Unrecognized file format.", type)
        return


# Count single nucleotide polymorphisms.
def get_pairwise_differences(first_sequence, second_sequence):
    return np.sum(x != y for x, y in zip(first_sequence, second_sequence))\
            + abs(len(first_sequence) - len(second_sequence))


# Calculate nCr.
def comb(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok//ktok
    else:
        return 0


# Calculate nucleotide diversity.
def get_nucleotide_diversity(sequences):
    pi_value = 0
    for first_sequence in sequences:
        for second_sequence in sequences:
            pi_value += get_pairwise_differences(first_sequence, second_sequence)
    return pi_value/(2*comb(len(sequences), 2))


# Record segregating sites.
def get_segregating_sites(sequences):
    sites = []
    for first_sequence in sequences:
        for second_sequence in sequences:
            upper_boundary = min(len(first_sequence), len(second_sequence))
            current_list = [index for index in range(upper_boundary) if
                            first_sequence[index] != second_sequence[index]]
            sites = sorted(np.unique(sites + current_list))
    return sites


# Calculate Tajima's D.
def get_tajimas_d(sequences):
    n = len(sequences)
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
    s = len(get_segregating_sites(sequences))
    k = get_nucleotide_diversity(sequences)
    avg_length = 0
    for sequence in sequences:
        avg_length += len(sequence)/n
    if (np.sqrt(e_1*s + e_2*s*(s - 1))) == 0:
        return float("nan")
    else:
        return (k - s/a_1)/((np.sqrt(e_1*s + e_2*s*(s - 1)))*avg_length)


# Parse sequence into pieces with fixed length.
def parse_into_pieces(sequences, window_size):
    pieces = []
    index_1 = 0
    for sequence in sequences:
        num = len(sequence)//window_size
        pieces = ([None]*len(sequences))
        for index_2 in range(num):
            new_piece = sequence[index_2*window_size:(index_2 + 1)*window_size]
            if new_piece:
                pieces[index_1].append(new_piece)
        if sequence[num*window_size:]:
            pieces[index_1].append(sequence[num*window_size:])
        index_1 += 1
    return pieces


# Calculate the Tajima's D for each piece
def analyze_pieces(parsed_sequences, window_size):
    positions = []
    scores = []
    for column in range(len(parsed_sequences[0])):
        current = []
        for row in range(len(parsed_sequences)):
            current.append(parsed_sequences[row][column])
        if not np.isnan(get_tajimas_d(current)):
            scores.append(get_tajimas_d(current))
            if len(current[0]) < window_size:
                positions.append(column*window_size + (len(current[0]) + 1)//2)
            else:
                positions.append(column*window_size + (window_size + 1)//2)
    return [positions, scores]


# The main function.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Balancing selection detector.')
    parser.add_argument(dest='file_name', help='Please enter the filename.')
    parser.add_argument(dest='window_size', help='Please enter the window size.', type=int)
    args = parser.parse_args()

    print("The input file is: " + str(args.file_name))
    if os.path.isfile("sequence_list.pkl"):
        pkl_file = open(f"{args.file_name}.pkl", "rb")
        Sequences = pickle.load(pkl_file)
    else:
        Sequences = read_sequences(args.file_name)
        with open(f"{args.file_name}.pkl", "wb") as p:
            pickle.dump(Sequences, p)

    print("The number of sample is: " + str(len(Sequences)))
    # print("The sequences are: ")
    # ppt.pprint(Sequences)

    Pi_value = get_nucleotide_diversity(Sequences)
    print("The nucleotide diversity is: " + str(Pi_value))
    Tajima_D = get_tajimas_d(Sequences)
    print("The Tajima's D is: " + str(Tajima_D))

    print("The input window size is: " + str(args.window_size))
    Parsed_sequences = parse_into_pieces(Sequences, args.window_size)
    # print("The parsed sequences are: ")
    # ppt.pprint(Parsed_sequences)

    [Bp_positions, Tajima_scores] = analyze_pieces(Parsed_sequences, args.window_size)
    plt.plot(Bp_positions, Tajima_scores, color='blue', linestyle='dashed', linewidth=1,
             marker='.', markerfacecolor='blue', markersize=5)
    plt.xlim(0, max(Bp_positions))
    plt.ylim(-1.5*max(Tajima_scores), 1.5*max(Tajima_scores))
    plt.xlabel("Position")
    plt.ylabel("Tajima's D")
    plt.title(f"{args.file_name} Balancing Selection Analysis")
    plt.savefig(f"{args.file_name} Tajima's D.png")
