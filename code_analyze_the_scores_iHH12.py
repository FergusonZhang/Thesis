# This program will determine the cutoff iHH12 value for the whole genome
# It will then plot iHH12 vs. Position figures for all eight chromosomes
# Finally, it will return the outliers with corresponding positions for each chromosome
import argparse
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from sklearn import preprocessing


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gene selector for iHH12.')
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    args = parser.parse_args()

    Data = pd.DataFrame()
    for i in range(1, 9):
        data = pd.read_csv(f'Results_iHH12/outfile_{i}.ihh12.out', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        data['d'] = preprocessing.scale(data['d'])
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
    with open(f'Results_iHH12/outliers_1.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_1, p)
        p.close()
    with open(f'Results_iHH12/outliers_2.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_2, p)
        p.close()
    with open(f'Results_iHH12/outliers_3.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_3, p)
        p.close()
    with open(f'Results_iHH12/outliers_4.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_4, p)
        p.close()
    with open(f'Results_iHH12/outliers_5.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_5, p)
        p.close()
    with open(f'Results_iHH12/outliers_6.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_6, p)
        p.close()
    with open(f'Results_iHH12/outliers_7.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_7, p)
        p.close()
    with open(f'Results_iHH12/outliers_8.ihh12.out.pkl', 'wb') as p:
        pickle.dump(Outliers_8, p)
        p.close()

    for j in range(1, 9):
        data = pd.read_csv(f'Results_iHH12/outfile_{j}.ihh12.out', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        Positions = data['b'].values.tolist()
        Scores = data['d'].values.tolist()
        infile = open(f'Results_iHH12/outliers_{j}.ihh12.out.pkl', 'rb')
        Outlier_positions = pickle.load(infile)
        infile.close()

        plt.figure(figsize=(40, 5))
        plt.plot(Positions, Scores, Color='blue', linewidth=0.5)
        # plt.plot(Outlier_positions, Outlier_scores, 'ro', markersize=2)
        plt.title('Balancing Selection Candidate Sites')
        plt.xlabel('Base Pair Position')
        plt.ylabel("iHH12")
        plt.savefig(f'Figures/Cgrand_scaffold_{j}_shapeit4.vcf_iHH12.png', dpi=500)
