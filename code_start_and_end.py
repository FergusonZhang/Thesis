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
        Genes = [gene - 4000 for gene in Genes]
        for gene in Genes:
            if gene <= Positions[0]:
                Genes.remove(gene)

        Gene_scores = np.zeros(4000)
        fragment_score = np.zeros(4000)
        for index, position in enumerate(Positions):
            if index < (len(Positions) - 4000):
                for gene in Genes:
                    if 0 <= (position - gene) <= 25:
                        point = index
                        while Positions[point] <= (gene + 4000):
                            location = int(Positions[point] - gene)
                            fragment_score[location] = Scores[point]
                            point += 1
                    zipped_lists = zip(Gene_scores, fragment_score)
                    Gene_scores = [x + y for (x, y) in zipped_lists]
                    fragment_score = np.zeros(4000)
        Gene_scores = [number/4000 for number in Gene_scores]

        x_values = list(range(-2000, 2000))
        plt.plot(x_values, Gene_scores, 'bo', markersize=1)
        plt.xlabel('Position')
        plt.ylabel("Tajima's D")
        plt.savefig(f'Data_extra/Chromosome_{i}_{args.name}')
