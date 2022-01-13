# This program will plot the iHH12 vs. Position figure as well as return outliers for each chromosome
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from sklearn import preprocessing


# The main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Outlier analyzer for iHH12.')
    parser.add_argument(dest='top', help='Please enter the top number.', type=int)
    args = parser.parse_args()

    Data = pd.DataFrame()
    for i in range(1, 9):
        data = pd.read_csv(f'Data_iHH12/Cgrand_scaffold_{i}_shapeit4.vcf.ihh12.out', sep='\t', header=None)
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
            Outliers_1.append([row['b'], row['d']])
        elif row['a'][9] == '2':
            Outliers_2.append([row['b'], row['d']])
        elif row['a'][9] == '3':
            Outliers_3.append([row['b'], row['d']])
        elif row['a'][9] == '4':
            Outliers_4.append([row['b'], row['d']])
        elif row['a'][9] == '5':
            Outliers_5.append([row['b'], row['d']])
        elif row['a'][9] == '6':
            Outliers_6.append([row['b'], row['d']])
        elif row['a'][9] == '7':
            Outliers_7.append([row['b'], row['d']])
        elif row['a'][9] == '8':
            Outliers_8.append([row['b'], row['d']])

    with open('Data_iHH12/Cgrand_scaffold_1_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_1, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_2_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_2, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_3_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_3, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_4_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_4, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_5_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_5, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_6_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_6, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_7_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_7, p)
        p.close()
    with open('Data_iHH12/Cgrand_scaffold_8_shapeit4.vcf_ihh12_candidates.pkl', 'wb') as p:
        pickle.dump(Outliers_8, p)
        p.close()

    for j in range(1, 9):
        data = pd.read_csv(f'Data_iHH12/Cgrand_scaffold_{j}_shapeit4.vcf.ihh12.out', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        data['d'] = preprocessing.scale(data['d'])
        Positions = data['b'].values.tolist()
        Scores = data['d'].values.tolist()
        infile = open(f'Data_iHH12/Cgrand_scaffold_{j}_shapeit4.vcf_ihh12_candidates.pkl', 'rb')
        Temp = pickle.load(infile)
        infile.close()
        Outlier_positions = [x[0] for x in Temp]
        Outlier_scores = [y[1] for y in Temp]

        # Plot the iHH12 vs. Position figure
        plt.figure(figsize=(25, 5))
        plt.plot(Positions, Scores, Color='blue', linewidth=0.5)
        plt.plot(Outlier_positions, Outlier_scores, 'ro', markersize=2)
        plt.title('Balancing Selection Candidate Sites')
        plt.xlabel('Position')
        plt.ylabel("iHH12")
        plt.savefig(f'Data_iHH12/Cgrand_scaffold_{j}_shapeit4.vcf_ihh12_figure.png', dpi=500)
