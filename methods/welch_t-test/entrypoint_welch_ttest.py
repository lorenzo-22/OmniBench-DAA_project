#!/usr/bin/env python3
"""
Welch T-Test Entrypoint for OmniBenchmark
"""

import argparse
import sys
import os

# Import the original welch_t-test script
sys.path.insert(0, os.path.dirname(__file__))

# Import functions from the original script
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def parse_arguments():
    parser = argparse.ArgumentParser(description='Welch t-test for OmniBenchmark')
    
    # OmniBenchmark standard arguments
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--name', type=str, required=True,
                        help='Dataset name')
    parser.add_argument('--data.matrix', dest='data_matrix', type=str, required=True,
                        help='Input data file')
    
    # Method parameters with defaults for this dataset
    parser.add_argument('--group1', default='R,S,T',
                        help='Comma-separated column names for group 1')
    parser.add_argument('--group2', default='U,V,W',
                        help='Comma-separated column names for group 2')
    parser.add_argument('--sheet', default=0,
                        help='Sheet name or index')
    parser.add_argument('--fdr', type=float, default=0.05,
                        help='FDR threshold')
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Output file path
    output_file = os.path.join(args.output_dir, f'{args.name}_results.csv')
    
    print(f"Running Welch t-test on {args.data_matrix}")
    print(f"Group 1: {args.group1}")
    print(f"Group 2: {args.group2}")
    print(f"Output: {output_file}")
    
    # Load data
    df = pd.read_excel(args.data_matrix, sheet_name=args.sheet)
    
    # Parse groups
    group1_cols = [c.strip() for c in args.group1.split(',')]
    group2_cols = [c.strip() for c in args.group2.split(',')]
    
    # Perform t-test (simplified version)
    results = []
    for idx, row in df.iterrows():
        g1_vals = row[group1_cols].dropna().values
        g2_vals = row[group2_cols].dropna().values
        
        if len(g1_vals) >= 2 and len(g2_vals) >= 2:
            t_stat, p_value = ttest_ind(g2_vals, g1_vals, equal_var=False)
            results.append({
                'Feature': idx,
                'P_Value': p_value,
                'T_Statistic': t_stat
            })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    if len(results_df) > 0:
        results_df['P_Adjusted'] = multipletests(results_df['P_Value'], method='fdr_bh')[1]
        results_df['Significant'] = results_df['P_Adjusted'] < args.fdr
    
    # Save results
    results_df.to_csv(output_file, index=False)
    
    print(f"✓ Results saved to: {output_file}")
    print(f"  Total features: {len(results_df)}")
    if len(results_df) > 0:
        print(f"  Significant: {results_df['Significant'].sum()}")


if __name__ == "__main__":
    main()
