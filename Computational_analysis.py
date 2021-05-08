# This file contains computational analysis methods
import argparse
import numpy as np
import warnings
import pprint as ppt
warnings.filterwarnings("ignore")


# Read sequences from different files
def read_sequences(file_name):
    file_type = file_name.split(".")[1]
    if file_type == "fas":
        sequences = []
        with open(file_name) as f:
            current_sequence = ""
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
        return
    else:
        print("Unrecognized file format!", type)
        return


# Count single nucleotide polymorphism between two sequences.
def pairwise_difference(first_sequence, second_sequence):
    return np.sum(x != y for x, y in zip(first_sequence, second_sequence))\
            + abs(len(first_sequence) - len(second_sequence))


# Calculate nCr without using the math package
def comb(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


# Calculate nucleotide diversity using pairwise difference.
def nucleotide_diversity(sequences):
    pi_value = 0
    for first_sequence in sequences:
        for second_sequence in sequences:
            pi_value += pairwise_difference(first_sequence, second_sequence)
    return pi_value/(2*comb(len(sequences), 2))


# Count the number of segregating site.
def segregating_site(sequences):
    result_list = []
    for first_sequence in sequences:
        for second_sequence in sequences:
            upper_boundary = min(len(first_sequence), len(second_sequence))
            current_list = [index for index in range(upper_boundary) if
                            first_sequence[index] != second_sequence[index]]
            result_list = sorted(np.unique(result_list + current_list))
    return len(result_list)


# Calculate Tajima's D using nucleotide diversity.
def tajima_test(sequences):
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
    s = segregating_site(sequences)
    k = nucleotide_diversity(sequences)
    avg_length = 0
    for sequence in sequences:
        avg_length += len(sequence)/n
    return (k - s/a_1)/((np.sqrt(e_1*s + e_2*s*(s - 1)))*avg_length)


# Parse sequence into pieces with fixed length
def parse_into_pieces(sequences, window_size):
    pieces = []
    for index_1, sequence in enumerate(sequences):
        num = len(sequence)//window_size
        pieces.append([]*num)
        for index_2 in range(num):
            new_piece = sequence[index_2*window_size:(index_2 + 1)*window_size]
            if new_piece:
                pieces[index_1].append(new_piece)
        if sequence[num*window_size:]:
            pieces[index_1].append(sequence[num * window_size:])
    return pieces


# The main function.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Balancing selection detector.')
    parser.add_argument(dest='file_name', help='Please enter the filename.')
    parser.add_argument(dest='Window_size', help='Please enter the window size.', type=int)
    args = parser.parse_args()

    print("The input file is: " + str(args.file_name))
    Sequences = read_sequences(args.file_name)

    print("The sequences are: ")
    ppt.pprint(Sequences)

    Pi_value = nucleotide_diversity(Sequences)
    print("The nucleotide diversity is: " + str(Pi_value))

    Tajima_D = tajima_test(Sequences)
    print("The Tajima's D is: " + str(Tajima_D))

    Parsed_sequences = parse_into_pieces(Sequences, args.Window_size)
    print("The parsed sequences are: ")
    ppt.pprint(Parsed_sequences)

    Bp_position = []
    Tajima_scores = []
    for i in range(len(Parsed_sequences[0])):
        current_column = []
        for j in range(len(Parsed_sequences)):
            current_column.append(Parsed_sequences[j][i])
        if len(current_column[0]) < args.Window_size:
            Bp_position.append(i*args.Window_size + (len(current_column[0]) + 1)//2)
        else:
            Bp_position.append(i*args.Window_size + (args.Window_size + 1)//2)
        Tajima_scores.append(tajima_test(current_column))
    print("The score vs. position is: ")
    ppt.pprint(Bp_position)
    ppt.pprint(Tajima_scores)
