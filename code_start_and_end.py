# An interesting extra work
import argparse
import gffutils
import matplotlib.pyplot as plt
import numpy as np
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Start or end analysis.')
    parser.add_argument(dest='name', help='Please enter start or end.')
    args = parser.parse_args()
    db = gffutils.create_db('Data_raw/Crubella_474_v1.1.gene.gff3', 'data.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)
    for i in range(1, 2):
        infile = open(f'Data_extra/Cgrand_scaffold_{i}_shapeit4.vcf_positions.pkl', 'rb')
        Positions = pickle.load(infile)
        infile.close()
        infile = open(f'Data_extra/Cgrand_scaffold_{i}_shapeit4.vcf_scores.pkl', 'rb')
        Scores = pickle.load(infile)
        infile.close()

        Genes = []
        if args.name == 'start':
            for mRNA in db.features_of_type('mRNA', order_by='start'):
                if mRNA.seqid == f'scaffold_{i}':
                    Genes.append(mRNA.start)
        elif args.name == 'end':
            for mRNA in db.features_of_type('mRNA', order_by='end'):
                if mRNA.seqid == f'scaffold_{i}':
                    Genes.append(mRNA.end)
        Gene_scores = np.zeros(4000)
        for index, position in enumerate(Positions):
            for gene in Genes:
                if (index > 2000) and (abs(position - gene) <= 10) and (index < len(Positions) - 2001):
                    fragment_position = Positions[(index - 2000):(index + 2001)]
                    fragment_score = Scores[(index - 2000):(index + 2001)]
                    for point, element in enumerate(fragment_position):
                        if (point < position - 2000) or (point > position + 2000):
                            fragment_score[point] = 0
                    zipped_lists = zip(Gene_scores, fragment_score)
                    Gene_scores = [x + y for (x, y) in zipped_lists]
        Gene_scores = [number/4000 for number in Gene_scores]
        x_values = list(range(1, len(Gene_scores) + 1))
        plt.plot(x_values, Gene_scores, 'bo', markersize=1)
        plt.xlabel('Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Chromosome_{i}_{args.name}', dpi=500)
