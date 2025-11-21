#!/usr/bin/env python3
"""
Welch T-Test Calculator
Performs Welch t-test on proteomics data from command line

Usage:
    python welch_ttest.py -i input.csv -o results.csv --group1 A1,A2,A3 --group2 B1,B2,B3
"""

import argparse
import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Perform Welch t-test on quantitative data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with CSV file
  python welch_ttest.py -i data.csv -o results.csv --group1 A1,A2,A3 --group2 B1,B2,B3
  
  # Excel file with sheet specification
  python welch_ttest.py -i data.xlsx -o results.csv --sheet "Intensity MaxQuant (6)" \\
                        --group1 "A1 (log),A2 (log),A3 (log)" \\
                        --group2 "B1 (log),B2 (log),B3 (log)"
  
  # With ID column and FDR correction
  python welch_ttest.py -i data.csv -o results.csv \\
                        --group1 Low_1,Low_2,Low_3 \\
                        --group2 High_1,High_2,High_3 \\
                        --id-col Protein_ID \\
                        --fdr 0.05
  
  # TSV file with custom separator
  python welch_ttest.py -i data.tsv -o results.tsv --sep "\\t" \\
                        --group1 A1,A2,A3 --group2 B1,B2,B3
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Input file path (CSV, TSV, or Excel)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path')
    parser.add_argument('--group1', required=True,
                        help='Comma-separated column names for group 1 (e.g., A1,A2,A3)')
    parser.add_argument('--group2', required=True,
                        help='Comma-separated column names for group 2 (e.g., B1,B2,B3)')
    
    # Optional arguments
    parser.add_argument('--id-col', default=None,
                        help='Column name for feature IDs (default: first column)')
    parser.add_argument('--sheet', default=0,
                        help='Sheet name or index for Excel files (default: 0)')
    parser.add_argument('--sep', default=',',
                        help='Separator for CSV/TSV files (default: ",")')
    parser.add_argument('--fdr', type=float, default=0.05,
                        help='FDR threshold for significance (default: 0.05)')
    parser.add_argument('--min-valid', type=int, default=2,
                        help='Minimum valid values required per group (default: 2)')
    parser.add_argument('--log-transform', action='store_true',
                        help='Apply log2 transformation to data')
    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed progress information')
    
    return parser.parse_args()

def load_data(args):
    """Load data from input file"""
    if args.verbose:
        print(f"Loading data from: {args.input}")
    
    # Determine file type and load
    if args.input.endswith('.xlsx') or args.input.endswith('.xls'):
        # Excel file
        try:
            sheet = int(args.sheet) if args.sheet.isdigit() else args.sheet
        except:
            sheet = args.sheet
        df = pd.read_excel(args.input, sheet_name=sheet)
    else:
        # CSV/TSV file
        df = pd.read_csv(args.input, sep=args.sep)
    
    if args.verbose:
        print(f"Loaded {len(df)} rows, {len(df.columns)} columns")
    
    return df

def parse_column_names(col_string):
    """Parse comma-separated column names"""
    return [col.strip() for col in col_string.split(',')]

def welch_ttest(df, args):
    """Perform Welch t-test"""
    
    # Parse column names
    group1_cols = parse_column_names(args.group1)
    group2_cols = parse_column_names(args.group2)
    
    if args.verbose:
        print(f"\nGroup 1 columns ({len(group1_cols)}): {group1_cols}")
        print(f"Group 2 columns ({len(group2_cols)}): {group2_cols}")
    
    # Verify columns exist
    missing_cols = []
    for col in group1_cols + group2_cols:
        if col not in df.columns:
            missing_cols.append(col)
    
    if missing_cols:
        print(f"\nERROR: The following columns were not found in the data:")
        for col in missing_cols:
            print(f"  - {col}")
        print("\nAvailable columns:")
        for i, col in enumerate(df.columns[:20], 1):
            print(f"  {i}. {col}")
        if len(df.columns) > 20:
            print(f"  ... and {len(df.columns) - 20} more")
        sys.exit(1)
    
    # Determine ID column
    if args.id_col:
        if args.id_col not in df.columns:
            print(f"ERROR: ID column '{args.id_col}' not found")
            sys.exit(1)
        id_col = args.id_col
    else:
        id_col = df.columns[0]
    
    if args.verbose:
        print(f"Using ID column: {id_col}")
    
    # Extract data
    working_df = df[[id_col] + group1_cols + group2_cols].copy()
    
    # Convert to numeric
    for col in group1_cols + group2_cols:
        working_df[col] = pd.to_numeric(working_df[col], errors='coerce')
    
    # Log transform if requested
    if args.log_transform:
        if args.verbose:
            print("Applying log2 transformation...")
        for col in group1_cols + group2_cols:
            working_df[col] = np.log2(working_df[col] + 1)
    
    # Filter based on minimum valid values
    group1_valid = working_df[group1_cols].notna().sum(axis=1)
    group2_valid = working_df[group2_cols].notna().sum(axis=1)
    
    valid_mask = (group1_valid >= args.min_valid) & (group2_valid >= args.min_valid)
    working_df = working_df[valid_mask].copy()
    
    if args.verbose:
        print(f"\nFiltering: keeping features with ≥{args.min_valid} valid values per group")
        print(f"Features after filtering: {len(working_df)}")
    
    if len(working_df) == 0:
        print("ERROR: No features passed filtering criteria")
        sys.exit(1)
    
    # Perform Welch t-test for each row
    if args.verbose:
        print("\nPerforming Welch t-tests...")
    
    results = []
    for idx, row in working_df.iterrows():
        feature_id = row[id_col]
        
        # Get values and ensure float type
        g1_vals = row[group1_cols].dropna().astype(float).values
        g2_vals = row[group2_cols].dropna().astype(float).values
        
        # Calculate statistics
        g1_mean = np.mean(g1_vals)
        g2_mean = np.mean(g2_vals)
        g1_std = np.std(g1_vals, ddof=1)
        g2_std = np.std(g2_vals, ddof=1)
        
        # Log2 fold change
        log2fc = g2_mean - g1_mean
        fold_change = 2 ** log2fc
        
        # Welch t-test
        if len(g1_vals) >= 2 and len(g2_vals) >= 2:
            t_stat, p_value = ttest_ind(g2_vals, g1_vals, equal_var=False)
        else:
            t_stat, p_value = np.nan, np.nan
        
        results.append({
            'ID': feature_id,
            'Group1_Mean': g1_mean,
            'Group1_SD': g1_std,
            'Group1_N': len(g1_vals),
            'Group2_Mean': g2_mean,
            'Group2_SD': g2_std,
            'Group2_N': len(g2_vals),
            'Log2FC': log2fc,
            'FoldChange': fold_change,
            'T_Statistic': t_stat,
            'P_Value': p_value
        })
    
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction
    if args.verbose:
        print("Applying FDR correction (Benjamini-Hochberg)...")
    
    valid_pvals = results_df['P_Value'].notna()
    results_df['P_Adjusted'] = np.nan
    
    if valid_pvals.sum() > 0:
        results_df.loc[valid_pvals, 'P_Adjusted'] = multipletests(
            results_df.loc[valid_pvals, 'P_Value'],
            method='fdr_bh'
        )[1]
    
    # Determine significance
    results_df['Significant'] = results_df['P_Adjusted'] < args.fdr
    
    # Sort by p-value
    results_df = results_df.sort_values('P_Value')
    
    return results_df

def save_results(results_df, args):
    """Save results to output file"""
    if args.verbose:
        print(f"\nSaving results to: {args.output}")
    
    # Determine output format
    if args.output.endswith('.xlsx'):
        results_df.to_excel(args.output, index=False)
    else:
        sep = '\t' if args.output.endswith('.tsv') else ','
        results_df.to_csv(args.output, sep=sep, index=False)
    
    if args.verbose:
        print(f"✓ Saved {len(results_df)} results")

def print_summary(results_df, args):
    """Print summary statistics"""
    print("\n" + "=" * 70)
    print("WELCH T-TEST RESULTS SUMMARY")
    print("=" * 70)
    
    total = len(results_df)
    significant = results_df['Significant'].sum()
    
    print(f"\nTotal features tested: {total}")
    print(f"Significant (FDR < {args.fdr}): {significant} ({100*significant/total:.1f}%)")
    
    if significant > 0:
        sig_df = results_df[results_df['Significant']]
        print(f"\nSignificant features statistics:")
        print(f"  Mean |Log2FC|: {np.abs(sig_df['Log2FC']).mean():.3f}")
        print(f"  Median p-value: {sig_df['P_Value'].median():.2e}")
        
        print(f"\nTop 10 most significant features:")
        print(sig_df[['ID', 'Log2FC', 'P_Value', 'P_Adjusted']].head(10).to_string(index=False))
    
    print("\n" + "=" * 70)

def main():
    """Main function"""
    args = parse_arguments()
    
    # Load data
    df = load_data(args)
    
    # Perform Welch t-test
    results_df = welch_ttest(df, args)
    
    # Save results
    save_results(results_df, args)
    
    # Print summary
    print_summary(results_df, args)
    
    print(f"\n✓ Analysis complete! Results saved to: {args.output}")

if __name__ == "__main__":
    main()
