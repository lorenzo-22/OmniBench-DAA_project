#!/usr/bin/env python

"""
Simple Welch t-test Omnibenchmark module.

Loads:
  --data         csv dataframe with proteins as rows and samples as columns. 

Outputs:
  A CSV-like matrix:
      header: protein, effect_size, p-value
      rows:   per-protein results
"""

import argparse
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import warnings


def load_dataset(data_file):
    data = (
    pd.read_csv(data_file, index_col=0)
    .select_dtypes(include=['number']) # Select only proteomics data
    .T #transform
)
    if data.ndim != 2:
        raise ValueError("Invalid data structure: data.matrix must be 2D")
    return data

def load_labels(labels_file):
    data = (
        pd.read_csv(labels_file, index_col=0)
    )
    labels = data.iloc[:, 0].to_numpy()
    return labels

# def get_labels(data):
#     # Make sure there is something in the index
#     if len(data.index) == 0:
#         raise ValueError("Data index is empty.")

#     labels = [0 if ".N" in name else 1 for name in data.index]

#     # Warn if there are only one class (case or control) in the data.
#     if len(set(labels)) < 2:
#         warnings.warn("Labels contain only one class.", UserWarning)
    
#     print(f'Found {labels.count(0)} controls and {labels.count(1)} cases in the data.')

#     return labels


def welch_ttest_df(data, labels):
    """
    data: Pandas dataframe with proteins as columns and samples as rows
    labels: 1D array of 0/1 group labels, length = n_samples = n_rows in data
    
    Returns:
        pandas DataFrame with t-statistics and p-values for each column
    """

    X = np.asarray(data)
    labels = np.asarray(labels)
    
    results = []

    for col_idx in range(data.shape[1]):
        col = X[:, col_idx]
        
        # Split into groups
        group0 = col[labels == 0]
        group1 = col[labels == 1]
        
        # Run Welch's t-test (unequal variance)
        t_stat, p_val = ttest_ind(group0, group1, equal_var=False)
        
        results.append({'protein': f'{data.columns[col_idx]}', 'effect_size': t_stat, 'p_value': p_val})

    df = pd.DataFrame(results)
    return df



def main():
    parser = argparse.ArgumentParser(description='Welch t-test benchmark runner')

    parser.add_argument('--data.matrix', type=str,
                        help='csv dataframe with proteins as rows and samples as columns.', required = True)
    parser.add_argument('--data.true_labels', type=str,
                        help='csv file with true labels.', required = True)
    parser.add_argument('--output_dir', type=str,
                        help='output directory to store data files.', 
                        required=True)
    # parser.add_argument('--name', type=str, help='name of the dataset', default='clustbench')
    # parser.add_argument('--method', type=str,
    #                     help='sklearn method',
    #                     required = True)

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    print('Loading data')
    data = load_dataset(getattr(args, 'data.matrix'))
    labels = load_labels(getattr(args, 'data.true_labels'))

    print('Running Welch t-test')
    results = welch_ttest_df(data, labels)

    results.to_csv(os.path.join(args.output_dir, "results_welch_ttest.csv"))
    print(f'Welch t-test results are stored in {os.path.join(args.output_dir, "results_welch_ttest.csv")}')

if __name__ == "__main__":
    main()