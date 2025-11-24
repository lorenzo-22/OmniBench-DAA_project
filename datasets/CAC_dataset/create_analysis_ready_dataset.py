import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required = True, help="Input dataset to be manipulated")

parser.parse_args()


def main():

    args = parser.parse_args()
    input_data_path = args.input

    input_dataframe = pd.read_csv(input_data_path).loc[:, ~pd.read_csv(input_data_path).columns.str.contains('^Unnamed')]
    names_ensembl_ids = input_dataframe[["Name","Database_ID"]]
    # ready to use expression matrix
    expression_matrix = input_dataframe.drop(["Name", "Database_ID"], axis = 1)
    expression_matrix.to_csv("input_data_ready_for_analysis.csv", index = True)
    # labels with name to be mapped for analysis
    sample_ids = expression_matrix.columns#colnames(expression_m_samples_only)[2:length(colnames(expression_m_samples_only))]
    labels_dict = {}
    for id in sample_ids:
        if "N" in id:
            labels_dict[id] = "N"
        else:
            labels_dict[id] = "T"
    sample_labels_df = pd.DataFrame(list(labels_dict.items()), columns=["Sample ID", "Group Label"])
    sample_labels_df.to_csv("samples_labels_groups.csv", index=True)
    
    # dataframe containing protein Names, Protein ensembl ID & indexes 
    metadata_input_dataset = input_dataframe[["Name","Database_ID"]]
    metadata_input_dataset.to_csv("metadata_input_dataset.csv", index=True)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {e}")
        # Optionally, exit with a non-zero code
        import sys
        sys.exit(1)


#expression_m <- read.csv("datasets/CAC_dataset/coad_protein_expression_matrix.csv")
#head(expression_m, 3)
#names_ensembl_ids <- data.frame(expression_m$Name, expression_m$Database_ID)
#colnames(names_ensembl_ids) <- c("Name", "Database_ID")

#"expression_m_samples_only <- expression_m[, !names(expression_m) %in% c("Name","Database_ID")]
#head(expression_m_samples_only, 3)
#colnames(expression_m_samples_only) <- sub("^X", "", colnames(expression_m_samples_only))
#sample_ids = colnames(expression_m_samples_only)[2:length(colnames(expression_m_samples_only))]