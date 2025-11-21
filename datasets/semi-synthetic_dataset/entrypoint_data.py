#!/usr/bin/env python3
"""
Data module for OmniBenchmark
Provides the semi-synthetic dataset
"""

import argparse
import shutil
import os


def main():
    parser = argparse.ArgumentParser(description='Provide semi-synthetic dataset')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--name', type=str, required=True,
                        help='Dataset name')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Source file
    source_file = os.path.join(os.path.dirname(__file__), 'semi-synthetic_dataset.xlsx')
    
    # Destination file
    dest_file = os.path.join(args.output_dir, f'{args.name}.xlsx')
    
    # Copy the file
    shutil.copy2(source_file, dest_file)
    
    print(f"✓ Dataset '{args.name}' copied to: {dest_file}")


if __name__ == "__main__":
    main()
