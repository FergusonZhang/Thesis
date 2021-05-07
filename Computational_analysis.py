# This file contains the calculation of nucleotide diversity.
import numpy as np
import warnings
from Parse_the_input import *
from pprint import pprint
warnings.filterwarnings("ignore")


def pairwise_difference(first_sequence, second_sequence):
    count = np.sum(nucleotide_A != nucleotide_B for nucleotide_A, nucleotide_B
                   in zip(first_sequence, second_sequence))\
            + abs(len(first_sequence) - len(second_sequence))
    return count


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


def nucleotide_diversity(sequences):
    pi_value = 0
    for first_sequence in sequences:
        for second_sequence in sequences:
            pi_value += pairwise_difference(first_sequence, second_sequence)
    pi_value = pi_value/(2*comb(len(sequences), 2))
    return pi_value


def get_segregating_sites(sequences):
    # TODO: find the number of segregating sites
    return 4


def tajima_test(sequences):
    n = len(sequences)
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
    s = get_segregating_sites(sequences)
    k = nucleotide_diversity(sequences)
    tajima_score = (k - s/a_1)/np.sqrt(e_1*s + e_2*s*(s - 1))
    return tajima_score


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genome reader.')
    parser.add_argument(dest='file_name', help='Please input the genome filename.')
    args = parser.parse_args()
    print("The input file is: " + str(args.file_name))
    Sequences = file_parser(args.file_name)
    print("The sequences are: ")
    pprint(Sequences)
    Pi_value = nucleotide_diversity(Sequences)
    print("The nucleotide diversity is: " + str(Pi_value))
    Tajima_D = tajima_test(Sequences)
    print("The Tajima's D is: " + str(Tajima_D))
