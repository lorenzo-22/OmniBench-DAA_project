### Datasets introduction

## Colon Adenocarcinoma Datasets
The dataset has been download from the CPTAC Data Portal using the cptac Python API. (NOTE: The script is not included in the repo for now).
The .csv file contains the abundance matrix of proteins (rows) across samples (columns). The abundance values have been obtained from MS-TMT based experimental technique and have been normalized and log2transformed by the University of Michigan Pipeline (umich).
The folder contains 2 versions of the Colon Adenocarcinoma dataset (COAD): 
- `abundance_matrix_coad_umich_with_NAN.csv` version contains NAN values; It contains 9457 proteins abundance values (includes proteins with NAN values)  measured across 197 samples (100 Normal, 97 Tumor)
- `abundance_matrix_coad_umich.csv` version has been filtered to exclude NAN values; It contains 4943 proteins abundance values measured across 197 samples (100 Normal, 97 Tumor)

Both the datasets have been filtered a priori to remove uniprot unreviewed protein entries. 
