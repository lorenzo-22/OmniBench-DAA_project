options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
install.packages("dplyr")
install.packages("tibble")

library(BiocManager)
library(limma)
library(dplyr)
library(tibble)

expression_m <- read.csv("../../datasets/CAC_dataset/coad_protein_expression_matrix.csv")
head(expression_m, 3)
names_ensembl_ids <- data.frame(expression_m$Name, expression_m$Database_ID)
colnames(names_ensembl_ids) <- c("Name", "Database_ID")

expression_m_samples_only <- expression_m[, !names(expression_m) %in% c("Name","Database_ID")]
head(expression_m_samples_only, 3)
colnames(expression_m_samples_only) <- sub("^X", "", colnames(expression_m_samples_only))
sample_ids = colnames(expression_m_samples_only)[2:length(colnames(expression_m_samples_only))]
labels = vector("list", length(sample_ids))
for (sid in 1:length(sample_ids))
{
    string_spl = strsplit(sample_ids[[sid]], "\\.")
    last_str = string_spl[[1]][length(string_spl[[1]])]
    if(last_str == "N")
    {
        labels[[sid]] <- last_str
    }
    else
    {
        labels[[sid]] <- "T"
    }
}

label_samples_df <- data.frame(sample_id = I(sample_ids), label = I(labels))
expr_matrix_limma_tmp <- expression_m[, !(colnames(expression_m) %in% c("Database_ID"))]
protein_names <- expression_m$Name
print(protein_names)
expr_matrix_limma_tmp1 <- expr_matrix_limma_tmp[, !(colnames(expr_matrix_limma_tmp) %in% c("Name"))]
expr_matrix_limma <- expr_matrix_limma_tmp1[, -1]
colnames(expr_matrix_limma) <- sub("^X", "", colnames(expr_matrix_limma))


# --- create group containing the explanotory variable for the DEA analysis
# -- the variables are labels for tumor and normal samples
labels_vec <- unlist(label_samples_df[,2])
groups <- factor(labels_vec)

# create the design matrix. In limma, the design matrix defines the design of the model and contains the information about the parameters to take into account for the linear model to be fit on the data later on
# here we specify a design where parameters to consider are related to the two groups (T vs N) of samples-without the intercept (meaning that we consider the means of the two groups independetly to be compared for DEA)
# we are saying: "we want to see if the gene is differentially expressed in one of the two groups. To do so, we want to compare for a each gene the mean across normal samples with the mean across tumor samples and calculate the difference"
design <- model.matrix(~0 + groups)
# fit the linear model to the expression matrix data of each protein/gene
fit <- lmFit(expr_matrix_limma, design) 

# Now that the model has estimated the groups means for each gene/protein in the dataset- we can proceed with the contrast matrix
# the contrast matrix is used when no intercept is defined in the design matrix meaning that the linear model is not calculating the differences across the two groups (no reference given)- but just estimating the two means
# to actually get information on these estimation, we need to define a rule that the model can use to compare the two groups.
# For example, in this simple case we are just interested in computing the difference in mean expression between tumor and normal samples
contrasts_matrix <- makeContrasts(groupsT - groupsN, levels = design)
# fit the contrast matrix
fit2 <- contrasts.fit(fit, contrasts_matrix)

# --- Empirical Bayes moderated t-test
# -- used to improve the accuracy of the variance estimates for each feature (gene/protein)
# -- this is usually done when the number of sample is low and the estimated variability for a feature can be the result of random fluctuations or instability due to small sample size
# -- the method compute a prior distribution of the variances of all the features and then shrink each individual variance towards the estimated prior distribtion
fit3 <- eBayes(fit2, robust=TRUE)
limma_results <- topTable(fit3, coef = 1, adjust.method="BH", number=Inf, sort.by = "none") 
limma_results <- limma_results %>%
  mutate(Name = protein_names) %>%   # map names based on original row order
  select(Name, everything()) 

head(limma_results,3)