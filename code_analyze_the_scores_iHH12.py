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

    Data = pd.DataFrame()
    for i in range(1, 9):
        data = pd.read_csv(f'Results_iHH12/outfile_{i}.ihh12.out', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        Data = pd.concat([Data, data])
    Sorted_data = Data.sort_values(by=['d'], ascending=False, ignore_index=True)
    Sorted_data = Sorted_data.head(args.top)

    Outliers_1 = []
    Outliers_2 = []
    Outliers_3 = []
    Outliers_4 = []
    Outliers_5 = []
    Outliers_6 = []
    Outliers_7 = []
    Outliers_8 = []
    for index, row in Sorted_data.iterrows():
        if row['a'][9] == '1':
            Outliers_1.append(row['b'])
        elif row['a'][9] == '2':
            Outliers_2.append(row['b'])
        elif row['a'][9] == '3':
            Outliers_3.append(row['b'])
        elif row['a'][9] == '4':
            Outliers_4.append(row['b'])
        elif row['a'][9] == '5':
            Outliers_5.append(row['b'])
        elif row['a'][9] == '6':
            Outliers_6.append(row['b'])
        elif row['a'][9] == '7':
            Outliers_7.append(row['b'])
        elif row['a'][9] == '8':
            Outliers_8.append(row['b'])
    print(Outliers_2)

    for j in range(1, 9):
        data = pd.read_csv(f'Results_iHH12/outfile_{j}.ihh12.out', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        Positions = data['b'].values.tolist()
