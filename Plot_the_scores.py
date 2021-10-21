# This program will take a list of Tajima's D and a list of position to create a figure
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle

# The main function.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figure plotting.')
    parser.add_argument(dest='position_file_name', help='Please enter the position file name.')
    parser.add_argument(dest='score_file_name', help='Please enter the score file name.')
    parser.add_argument(dest='expected_size', help='Please enter the expected size.', type=int)
    args = parser.parse_args()

    infile = open(args.position_file_name, 'rb')
    Parsed_positions = pickle.load(infile)
    infile = open(args.score_file_name, 'rb')
    Tajimas_ds = pickle.load(infile)
    infile.close()

    index = np.argsort(Tajimas_ds)
    index = index[::-1]
    Sorted_positions = np.array(Parsed_positions)[index]
    Sorted_scores = np.array(Tajimas_ds)[index]
    Sorted_positions = Sorted_positions[0:args.expected_size]
    Sorted_scores = Sorted_scores[0:args.expected_size]

    plt.figure(figsize=(40, 5))
    plt.plot(Parsed_positions, Tajimas_ds, Color='blue', linewidth=0.5)
    plt.plot(Sorted_positions, Sorted_scores, 'ro', markersize=2)
    plt.title('Balancing Selection Candidate sites')
    plt.xlabel('Base Pair Position')
    plt.ylabel("Tajima's D")
    plt.savefig(f'{args.score_file_name}_figure.png', dpi=500)
