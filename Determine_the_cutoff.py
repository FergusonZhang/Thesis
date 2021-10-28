# This program will return the cut off value given a window size and an expected number
import argparse
import pickle

# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Cut off value determination.')
    parser.add_argument(dest='window_size', help='Please enter the parsing window size.', type=int)
    parser.add_argument(dest='expected_number', help='Please enter the expected number.', type=int)
    args = parser.parse_args()

    Tajimas_ds = []
    for i in range(1, 9):
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_scores.pkl', 'rb')
        Scores = pickle.load(infile)
        infile.close()
        Tajimas_ds.append(Scores)
    Sorted_scores = sorted(Tajimas_ds, reverse=True)
    print('The cut off value is: ' + Sorted_scores[args.expected_number])
