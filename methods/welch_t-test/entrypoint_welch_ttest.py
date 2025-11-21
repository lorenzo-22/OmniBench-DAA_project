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
    parser.add_argument('--group1', default='17,18,19',
                        help='Comma-separated column names or indices for group 1')
    parser.add_argument('--group2', default='20,21,22',
                        help='Comma-separated column names or indices for group 2')
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
    print(f"Group 1 columns: {args.group1}")
    print(f"Group 2 columns: {args.group2}")
    
    # Load data
    df = pd.read_excel(args.data_matrix, sheet_name=args.sheet)
    
    print(f"Data shape: {df.shape}")
    print(f"Column names: {df.columns.tolist()}")
    
    # Parse group column specifications
    # Support both column names and column indices (e.g., "17,18,19" for R,S,T)
    group1_spec = [c.strip() for c in args.group1.split(',')]
    group2_spec = [c.strip() for c in args.group2.split(',')]
    
    # Convert to column names or indices
    group1_cols = []
    group2_cols = []
    
    for spec in group1_spec:
        if spec.isdigit():
            group1_cols.append(df.columns[int(spec)])
        else:
            group1_cols.append(spec)
    
    for spec in group2_spec:
        if spec.isdigit():
            group2_cols.append(df.columns[int(spec)])
        else:
            group2_cols.append(spec)
    
    print(f"Group 1 using columns: {group1_cols}")
    print(f"Group 2 using columns: {group2_cols}")
    
    # Perform t-test for each row (feature/protein)
    results = []
    for idx in range(len(df)):
        row = df.iloc[idx]
        
        # Get values for each group
        try:
            g1_vals = row[group1_cols].dropna().astype(float).values
            g2_vals = row[group2_cols].dropna().astype(float).values
        except KeyError as e:
            print(f"Warning: Could not find columns for row {idx}: {e}")
            continue
        
        if len(g1_vals) >= 2 and len(g2_vals) >= 2:
            t_stat, p_value = ttest_ind(g2_vals, g1_vals, equal_var=False)
            
            # Get feature ID (first column)
            feature_id = df.iloc[idx, 0] if df.shape[1] > 0 else idx
            
            results.append({
                'Feature': feature_id,
                'Group1_Mean': np.mean(g1_vals),
                'Group2_Mean': np.mean(g2_vals),
                'Log2FC': np.log2(np.mean(g2_vals) + 1) - np.log2(np.mean(g1_vals) + 1),
                'T_Statistic': t_stat,
                'P_Value': p_value
            })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    if len(results_df) > 0:
        results_df['P_Adjusted'] = multipletests(results_df['P_Value'], method='fdr_bh')[1]
        results_df['Significant'] = results_df['P_Adjusted'] < args.fdr
    
    # Save results
    results_df.to_csv(output_file, index=False)
    
    print(f"✓ Results saved to: {output_file}")
    print(f"  Total features tested: {len(results_df)}")
    if len(results_df) > 0:
        print(f"  Significant (FDR < {args.fdr}): {results_df['Significant'].sum()}")


if __name__ == "__main__":
    main()
