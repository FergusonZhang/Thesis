# This program will plot the Tajima's D vs. Position figure as well as return outliers for each chromosome
import argparse
import matplotlib.pyplot as plt
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Outlier analyzer for Tajima's D.")
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    args = parser.parse_args()

    Scores = []
    for i in range(1, 9):
        infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_scores.pkl', 'rb')
        Scores.extend(pickle.load(infile))
        infile.close()
    Sorted_scores = sorted(Scores, reverse=True)
    Cutoff = Sorted_scores[args.top]

    for i in range(1, 9):
        Outlier_positions = []
        Outlier_scores = []
        infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_positions.pkl', 'rb')
        Parsed_positions = pickle.load(infile)
        infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_scores.pkl', 'rb')
        Tajimas_ds = pickle.load(infile)
        infile.close()
        for j in range(len(Parsed_positions)):
            if Tajimas_ds[j] > Cutoff:
                Outlier_positions.append(Parsed_positions[j])
                Outlier_scores.append(Tajimas_ds[j])
        with open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_candidates.pkl', 'wb') as p:
            pickle.dump(Outlier_positions, p)
        p.close()

        # Plot the Tajima's D vs. Position figure
        plt.figure(figsize=(25, 5))
        plt.plot(Parsed_positions, Tajimas_ds, Color='blue', linewidth=0.5)
        plt.plot(Outlier_positions, Outlier_scores, 'ro', markersize=2)
        plt.title('Balancing Selection Candidate Sites')
        plt.xlabel('Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_figure.png', dpi=500)
