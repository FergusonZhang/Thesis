# This program will take the lists of Tajima's D and position to create a figure.
import argparse
import matplotlib.pyplot as plt
import pickle

# The main function.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tajima figure plotting.')
    parser.add_argument(dest='positions_file_name', help='Please enter the base pair position filename.')
    parser.add_argument(dest='scores_file_name', help="Please enter the Tajima's D file name.")
    args = parser.parse_args()
    print('The input position file is: ' + str(args.positions_file_name))
    print('The input score file is: ' + str(args.scores_file_name))

    infile = open(args.positions_file_name, 'rb')
    Parsed_positions = pickle.load(infile)
    infile = open(args.scores_file_name, 'rb')
    Tajimas_ds = pickle.load(infile)
    infile.close()

    plt.plot(Parsed_positions, Tajimas_ds, Color='blue', linewidth=0.7)
    plt.xlabel('Base Pair Position')
    plt.ylabel("Tajima's D")
    plt.savefig(f'{args.file_name}_{args.window_size}_figure.png', dpi=200)
