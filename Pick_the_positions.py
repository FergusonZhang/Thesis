# This program will take a list of Tajima's D and a list of position and return BS candidates
import argparse
import numpy as np
import pickle

# The main function.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Peaks picking.')
    parser.add_argument(dest='position_file_name', help='Please enter the position file name.')
    parser.add_argument(dest='score_file_name', help='Please enter the score file name.')
    parser.add_argument(dest='expected_size', help='Please enter the expected size.', type=int)
    args = parser.parse_args()
    print('The input base pair position file is: ' + str(args.position_file_name))
    print("The input Tajima's D file is: " + str(args.score_file_name))

    infile = open(args.position_file_name, 'rb')
    Parsed_positions = pickle.load(infile)
    infile = open(args.score_file_name, 'rb')
    Tajimas_ds = pickle.load(infile)
    infile.close()

    index = np.argsort(Tajimas_ds)
    index = index[::-1]
    Sorted_positions = np.array(Parsed_positions)[index]
    Sorted_scores = np.array(Tajimas_ds)[index]
    print(Sorted_positions[0:args.expected_size])
    print(Sorted_scores[0:args.expected_size])
