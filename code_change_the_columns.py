# This code is designed to refine plink maps
import numpy as np
import pandas as pd

if __name__ == '__main__':
    for i in range(1, 9):
        data = pd.read_csv(f'Data_plink/scaffold_{i}_plink.map', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']

        # Calculate genetic positions assuming a uniform recombination possibility
        data['a'] = [np.int64(i)] * data.shape[0]
        if i == 1:
            constant = 1/19624005
        elif i == 2:
            constant = 1/14106445
        elif i == 3:
            constant = 1/15044833
        elif i == 4:
            constant = 1/15045541
        elif i == 5:
            constant = 1/13734110
        elif i == 6:
            constant = 1/16617047
        elif i == 7:
            constant = 1/17328217
        else:
            constant = 1/13362172
        data['c'] = constant*data['d']
        data.to_csv(f'Data_plink/modified_scaffold_{i}_plink.map', header=None, index=None, sep='\t')
