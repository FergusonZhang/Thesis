# This program will determine the cutoff iHH12 value for the whole genome
# It will then plot iHH12 vs. Position figures for all eight chromosomes
# Finally, it will return the outliers with corresponding positions for each chromosome
import argparse
import matplotlib.pyplot as plt
import pickle
import pandas as pd


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gene selector for iHH12.')
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    args = parser.parse_args()

    Data = pd.DataFrame(['a', 'b', 'c', 'd'])
    for i in range(1, 9):
        data = pd.read_csv(f'Results_iHH12/outfile_{i}.ihh12', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        Data = pd.concat([Data, data])
    Data.sort_values(by=['d'], axis=0, ascending=False)
    Sorted_data = Data.head(args.top)
    print(Sorted_data)
