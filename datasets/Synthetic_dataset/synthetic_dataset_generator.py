#!/usr/bin/env python3
import numpy as np
import pandas as pd

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
    true_labels = np.zeros(proteins_n)
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
    # Generate the dataframe
    df = generate_data()
    
    # Print summary to console
    print(f"\nData Generated Successfully!")
    print(f"Shape: {df.shape}")
    print("\nFirst 5 rows:")
    print(df.head())
    
    # Save to CSV
    output_filename = "simulated_proteomics_data.csv"
    df.to_csv(output_filename)
    print(f"\nSaved full dataset to '{output_filename}'")
