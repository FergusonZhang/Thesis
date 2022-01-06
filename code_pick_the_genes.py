# This program will take candidate sites and return candidate genes
import argparse
import gffutils
import pickle


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gene selector.')
    parser.add_argument(dest='task_name', help='Please enter the task name.')
    args = parser.parse_args()

    db = gffutils.create_db('Data_raw/Crubella_474_v1.1.gene.gff3', 'data.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)
    for i in range(1, 9):
        if args.task_name == "Tajima":
            infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_candidates.pkl', 'rb')
            Candidates = pickle.load(infile)
        else:
            infile = open(f'Data_iHH12/Cgrand_scaffold_{i}_shapeit4.vcf_ihh12_candidates.pkl', 'rb')
            Temp = pickle.load(infile)
            Candidates = [x[0] for x in Temp]
        infile.close()
        print('The candidates are: ')
        for gene in db.features_of_type('gene', order_by='start'):
            if gene.seqid == f'scaffold_{i}':
                for candidate in Candidates:
                    if gene.start <= candidate <= gene.end:
                        print(gene.featuretype, gene.start, gene.end, gene.id)
