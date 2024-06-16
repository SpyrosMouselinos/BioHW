"""
4. (4 pts) Write a program enrichment_test.py that, for a CSV file in the previously described format and a list of
 identifiers (rows), conducts a Fisher's test for each column (with Bonferroni correction).
 Use the program to check whether groups of proteins labeled as A and B in the fragments file differ significantly
 in terms of found domains or enrichment in GO functions of the genes encoding them.
"""
import json
import pandas as pd
import scipy.stats as stats



def fisher_test(matrix, groupA, groupB):
    """
    Perform Fisher's Test for each column in the matrix.
    """
    if isinstance(groupA, list):
        pass
    else:
        groupA = list(groupA)

    if isinstance(groupB, list):
        pass
    else:
        groupB = list(groupB)

    p_values = []
    for column in matrix.columns:
        # Create the contingency table
        groupA_present = matrix.loc[groupA, column].sum()
        groupA_absent = len(groupA) - groupA_present
        groupB_present = matrix.loc[groupB, column].sum()
        groupB_absent = len(groupB) - groupB_present

        contingency_table = [[groupA_present, groupA_absent], [groupB_present, groupB_absent]]

        # Perform Fisher's Exact Test
        _, p_value = stats.fisher_exact(contingency_table)
        p_values.append(p_value)

    # Apply Bonferroni correction
    bonferroni_corrected_p_values = [min(p * len(p_values), 1.0) for p in p_values]

    results = pd.DataFrame({
        'Feature': matrix.columns,
        'P-Value': p_values,
        'Bonferroni Corrected P-Value': bonferroni_corrected_p_values
    })

    return results

def main(input_csv,
         groupA,
         groupB,
         output_csv):
    # Load the matrix
    matrix = pd.read_csv(input_csv, index_col=0)

    # Perform Fisher's Exact Test with Bonferroni correction
    results = fisher_test(matrix, groupA, groupB)

    # Save the results to a CSV file
    results.to_csv(output_csv, index=False)
    print(f"Enrichment test results have been written to {output_csv}")


def extract_identifiers(group_to_gene_dict):
    with open(group_to_gene_dict, 'r') as fin:
        data = json.load(fin)
        groupA_identifiers = []
        groupB_identifiers = []

        for query_name, seqs in data.items():
            if query_name.startswith("groupA"):
                for ii in seqs:
                    groupA_identifiers.append(ii)
            elif query_name.startswith("groupB"):
                for jj in seqs:
                    groupB_identifiers.append(jj)

    return set(groupA_identifiers), set(groupB_identifiers)


if __name__ == "__main__":
    association_file = '../data/associations.json'
    groupA_identifiers, groupB_identifiers = extract_identifiers(association_file)


    # Let's begin with the PFAM DOMAINS from question 2 #
    input_csv = "../data/pfam_domains.csv"
    output_csv = "../data/enrichment_results_pfam.csv"
    main(input_csv, groupA_identifiers, groupB_identifiers, output_csv)

    # Then test GO_ANNOTATIONS from question 2 #
    input_csv = "../data/go_annotations.csv"
    output_csv = "../data/enrichment_results_go.csv"
    main(input_csv, groupA_identifiers, groupB_identifiers, output_csv)
