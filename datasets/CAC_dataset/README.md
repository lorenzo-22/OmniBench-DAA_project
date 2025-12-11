### Datasets introduction

## Colon Adenocarcinoma Datasets
The dataset has been download from the CPTAC Data Portal using the cptac Python API. (NOTE: The script is not included in the repo for now).
The .csv file contains the abundance matrix of proteins (rows) across samples (columns). The abundance values have been obtained from MS-TMT based experimental technique and have been normalized and log2transformed by the University of Michigan Pipeline (umich).
The folder contains 2 versions of the Colon Adenocarcinoma dataset (COAD): 
- `coad_protein_expression_matrix_with_NAN.csv` version contains NAN values; It contains 9457 proteins abundance values (includes proteins with NAN values)  measured across 197 samples (100 Normal, 97 Tumor)
- `coad_protein_expression_matrix.csv` version has been filtered to exclude NAN values; It contains 4943 proteins abundance values measured across 197 samples (100 Normal, 97 Tumor)

Both the datasets have been filtered a priori to remove uniprot unreviewed protein entries. 

The folder contains the script `create_analysis_ready_dataset.py`used to manipulate the starting dataset to produce a ready-to-use dataset for the statistical analysis. The new dataset only contains expression values over all the samples (the Names and Database_ID columns have been removed). The script also produce two additional dataframe:
- the `samples_labels_groups.csv`containing for each sample id the lable group of belonging (Tumor and Normal)
- the `metadata_input_dataset.csv` containing for each protein the corresponding ensembl protein ID (Database ID)