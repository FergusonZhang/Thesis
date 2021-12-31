# This code is designed to refine the plink map generated by the vcf package
import pandas as pd
import numpy as np

if __name__ == '__main__':
    column_names = ['a', 'b', 'c', 'd']
    all_data = pd.DataFrame(columns=column_names)
    for i in range(1, 9):
        data = pd.read_csv(f'Data_plink/scaffold_{i}_plink.map', sep='\t', header=None)
        data.columns = ['a', 'b', 'c', 'd']
        data['a'] = [np.int64(i)] * data.shape[0]
        constant = 0.000000001
        data['c'] = constant*data['d']
        frames = [all_data, data]
        all_data = pd.concat(frames)
    all_data.to_csv('Data_plink/modified_scaffold_plink.map', header=None, index=None, sep='\t')