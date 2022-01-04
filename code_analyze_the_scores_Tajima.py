# This program will determine the cutoff Tajima's D value for the whole genome
# It will then plot the figure of Tajima's D vs. Position for all eight chromosomes
# Finally, it will return the outliers with corresponding positions for each chromosome
import argparse
import matplotlib.pyplot as plt
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Outlier analyzer.')
    parser.add_argument(dest='window_size', help='Please enter the parsing window size.', type=int)
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    args = parser.parse_args()

    # Calculate the cutoff value based on the whole genome
    Scores = []
    for i in range(1, 9):
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_scores.pkl', 'rb')
        Scores.extend(pickle.load(infile))
        infile.close()
    Sorted_scores = sorted(Scores, reverse=True)
    Cutoff = Sorted_scores[args.top]

    # Pick outliers and corresponding positions for each chromosome
    for i in range(1, 9):
        Outlier_positions = []
        Outlier_scores = []
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_positions.pkl', 'rb')
        Parsed_positions = pickle.load(infile)
        infile = open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_scores.pkl', 'rb')
        Tajimas_ds = pickle.load(infile)
        infile.close()
        for j in range(len(Parsed_positions)):
            if Tajimas_ds[j] > Cutoff:
                Outlier_positions.append(Parsed_positions[j])
                Outlier_scores.append(Tajimas_ds[j])
        with open(f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_{args.top}_candidates.pkl', 'wb') as p:
            pickle.dump(Outlier_positions, p)
        p.close()

        # Plot the figure of Tajima's D vs. Position
        plt.figure(figsize=(40, 5))
        plt.plot(Parsed_positions, Tajimas_ds, Color='blue', linewidth=0.5)
        plt.plot(Outlier_positions, Outlier_scores, 'ro', markersize=2)
        plt.title('Balancing Selection Candidate Sites')
        plt.xlabel('Base Pair Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Figures/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_figure.png', dpi=500)
