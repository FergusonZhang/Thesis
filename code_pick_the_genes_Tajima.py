# This program will take candidate sites and return candidate genes
import gffutils
import pickle


# The main function
if __name__ == '__main__':
    db = gffutils.create_db('Data_raw/Crubella_474_v1.1.gene.gff3', 'data.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)
    for i in range(1, 9):
        infile = open(f'Data_Tajima/Cgrand_scaffold_{i}_shapeit4.vcf_candidates.pkl', 'rb')
        Candidates = pickle.load(infile)
        infile.close()

        for gene in db.features_of_type('gene', order_by='start'):
            if gene.seqid == f'scaffold_{i}':
                for candidate in Candidates:
                    if gene.start <= candidate <= gene.end:
                        print(gene.featuretype, gene.start, gene.end, gene.id)
