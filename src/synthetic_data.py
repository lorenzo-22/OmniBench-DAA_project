#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse

def generate_data():
    """Simulates proteomics data with a signal in the first 200 proteins."""
    print("Initializing simulation...")
    
    # 1. Setup parameters
    np.random.seed(42)
    proteins_n = 2000
    nA, nB = 6, 6
    mu, sigma = 20, 0.5
    
    # 2. Generate random normal distribution data
    A = np.random.normal(mu, sigma, size=(proteins_n, nA))
    B = np.random.normal(mu, sigma, size=(proteins_n, nB))
    
    # 3. Inject signal (+1 log2 effect) into first 200 proteins of Group B
    B[:200] += 1.0
    
    # 4. Create label tracker (optional, but good for validation)
    true_labels = np.zeros(proteins_n, dtype=int)
    true_labels[:200] = 1
    
    # 5. formatting into DataFrame
    counts = np.hstack([A, B])
    samples = [f"A{i}" for i in range(nA)] + [f"B{i}" for i in range(nB)]
    proteins = [f"Prot{i}" for i in range(proteins_n)]
    
    df = pd.DataFrame(counts, index=proteins, columns=samples)
    
    # Add the truth label as a column for reference (optional)
    df['is_differentially_expressed'] = true_labels
    
    return df

if __name__ == "__main__":
    # Setup Argument Parser for OmniBenchmark
    parser = argparse.ArgumentParser(description="Generate synthetic proteomics data")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory")
    parser.add_argument("--name", type=str, required=True,
                        help="Dataset name")
    
    args = parser.parse_args()

    # Generate the dataframe
    df = generate_data()
    
    # Create output directory if needed
    import os
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Construct output file paths
    output_file = os.path.join(args.output_dir, f"{args.name}.synthetic_dataset.csv")
    labels_file = os.path.join(args.output_dir, f"{args.name}.true_labels.csv")
    
    # Print summary to console
    print(f"\nData Generated Successfully!")
    print(f"Shape: {df.shape}")
    print("\nFirst 5 rows:")
    print(df.head())
    
    # Save main dataset to CSV (without labels column)
    df_data = df.drop(columns=['is_differentially_expressed'])
    df_data.to_csv(output_file)
    print(f"\nSaved dataset to '{output_file}'")
    
    # Save labels separately
    labels_df = df[['is_differentially_expressed']]
    labels_df.to_csv(labels_file)
    print(f"Saved true labels to '{labels_file}'")
