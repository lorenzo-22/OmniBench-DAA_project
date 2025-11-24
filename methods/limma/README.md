# LIMMA FOR DIFFERENTIAL ABUNDANCE ANALYSIS
## Limma
The objective of limma analysis is to detect statistically significant abundant proteins in two different condition groups. In this case, it will estimate Tumor and Normal group mean abundance for each protein and test whether the group means are significantly different. To this end, limma implements eBayes moderated t-statistics which 1) estimate the variance for each gene; 2) the common variance trend across all genes; 3) Shrinks each gene's variance estimate toward the common trend; 4) compute the statistical significance of the estimated differences.
The method takes as input the protein abundance matrix and returns as output a data table containing the Differential Abundant Proteins, each associated with log2FoldChange, p-values (adjusted) and other score metrics. 
## Requirments
To run limma locally, a conda env with R has been created
```bash
conda create --name myenv r-base=4.3
```
```bash
conda activate myenv r-base=4.3
``` 
```python
Rscript limma_method.r
``` 

