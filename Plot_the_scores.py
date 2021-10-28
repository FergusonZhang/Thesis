# This program will create Tajima's D vs. Positions figures
import argparse
import matplotlib.pyplot as plt
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figure plotting.')
    parser.add_argument(dest='window_size', help='Please enter the parsing window size.', type=int)
    parser.add_argument(dest='cutoff', help='Please enter the cut off value.', type=float)
    args = parser.parse_args()

    for i in range(1, 9):
        Outlier_positions = []
        Outlier_scores = []
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_positions.pkl', 'rb')
        Parsed_positions = pickle.load(infile)
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_scores.pkl', 'rb')
        Tajimas_ds = pickle.load(infile)
        infile.close()

        for j in range(len(Parsed_positions) + 1):
            if Tajimas_ds[j] >= args.cuttoff:
                Outlier_positions.append(Parsed_positions[j])
                Outlier_scores.append(Tajimas_ds[j])

        plt.figure(figsize=(40, 5))
        plt.plot(Parsed_positions, Tajimas_ds, Color='blue', linewidth=0.5)
        plt.plot(Outlier_positions, Outlier_scores, 'ro', markersize=2)
        plt.title('Balancing Selection Candidate Sites')
        plt.xlabel('Base Pair Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Figures/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_figure.png', dpi=500)
