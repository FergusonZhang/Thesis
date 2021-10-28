# This program will return the outlier scores with corresponding positions
import argparse
import pickle

# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument(dest='window_size', help='Please enter the parsing window size.', type=int)
    parser.add_argument(dest='cutoff', help='Please enter the cutoff value.', type=float)
    args = parser.parse_args()

    for i in range(1, 9):
        Outlier_positions = []
        Outlier_scores = []
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_positions.pkl', 'rb')
        Parsed_positions = pickle.load(infile)
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_scores.pkl', 'rb')
        Tajimas_ds = pickle.load(infile)
        infile.close()

        for j in range(len(Parsed_positions)):
            if Tajimas_ds[j] >= args.cutoff:
                Outlier_positions.append(Parsed_positions[j])
        print(f'The outliers on chromosome {i} are: ' + str(Outlier_positions))
