# This program will take candidate positions and return corresponding candidate genes
import argparse
import gffutils
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Gene selector for Tajima's D.")
    parser.add_argument(dest='window_size', help='Please enter the parsing window size.', type=int)
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    args = parser.parse_args()

    db = gffutils.create_db('Data_raw/Crubella_474_v1.1.gene.gff3', 'data.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)
    for i in range(1, 9):
        infile = open(
            f'Data_pkl/Cgrand_scaffold_{i}_shapeit4.vcf_{args.window_size}_{args.top}_candidates.pkl', 'rb')
        Candidates = pickle.load(infile)
        infile.close()
        print(f'The candidate genes on chromosome {i} are:')

        # Compare outliers with the reference
        for gene in db.features_of_type('gene', order_by='start'):
            if gene.seqid == f'scaffold_{i}':
                for candidate in Candidates:
                    if gene.start <= candidate <= gene.end:
                        print(gene.featuretype, gene.start, gene.end, gene.id)
        print('\n')
