# This program will take a list of candidate position on a chromosome and return candidates
import argparse
import gffutils
import pickle

# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gene selector.')
    parser.add_argument(dest='window_size', help='Please enter the parsing window size.', type=int)
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    parser.add_argument(dest='scaffold', help='Please enter the scaffold number.', type=int)
    args = parser.parse_args()

    infile = open(f'Data_pkl/Cgrand_scaffold_{args.scaffold}_{args.window_size}_{args.top}_candidates.pkl', 'rb')
    Candidates = pickle.load(infile)
    infile.close()

    db = gffutils.create_db('Data_raw/Crubella_474_v1.1.gene.gff3', 'data.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)
    Genes = []
    for gene in db.features_of_type('gene', order_by='start'):
        if gene.seqid == f'scaffold_{args.scaffold}':
            Genes.append(gene)

    for gene in Genes:
        for fragment in db.children(gene, order_by='start'):
            for candidate in Candidates:
                if fragment.start <= candidate <= fragment.end:
                    print(fragment.featuretype, fragment.id)
