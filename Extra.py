# An interesting extra work
import argparse
import gffutils
import matplotlib.pyplot as plt
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='?')
    args = parser.parse_args()
    db = gffutils.create_db('Data_raw/Crubella_474_v1.1.gene.gff3', 'data.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)
    for i in range(1, 2):
        infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_positions.pkl', 'rb')
        Positions = pickle.load(infile)
        infile.close()
        infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_scores.pkl', 'rb')
        Scores = pickle.load(infile)
        infile.close()

        Starts = []
        Ends = []
        for mRNA in db.features_of_type('mRNA', order_by='start'):
            if mRNA.seqid == f'scaffold_{i}':
                Starts.append(mRNA.start)
                Ends.append(mRNA.end)
        print(len(Starts))
        print(len(Ends))

        Start_scores = []
        count = 0
        for index, position in enumerate(Positions):
            for start in Starts:
                if (index > 20) and (abs(position - start) <= 100) and (index < len(Positions) - 21):
                    fragment = Scores[(index - 20):(index + 21)]
                    zipped_lists = zip(Start_scores, fragment)
                    Start_scores = [x + y for (x, y) in zipped_lists]
                    count += 1
        Start_scores = [number/count for number in Start_scores]
        print(len(Start_scores))

        End_scores = []
        count = 0
        for index, position in enumerate(Positions):
            for end in Ends:
                if (index > 20) and (abs(position - end) <= 100) and (end < len(Positions) - 21):
                    fragment = Scores[(index - 20):(index + 20)]
                    zipped_lists = zip(End_scores, fragment)
                    End_scores = [x + y for (x, y) in zipped_lists]
                    count += 1
        End_scores = [number/count for number in End_scores]
        print(len(End_scores))

        x_values = list(range(1, len(Start_scores) + 1))
        plt.plot(x_values, Start_scores, 'bo', markersize=0.5)
        plt.xlabel('Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Chromosome_{i}_start', dpi=500)

        x_values = list(range(1, len(End_scores) + 1))
        plt.plot(x_values, End_scores, 'bo', markersize=0.5)
        plt.xlabel('Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Chromosome_{i}_end', dpi=500)
