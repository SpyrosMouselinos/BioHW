"""
Check whether groups of proteins labeled as A and B in the fragments file differ significantly
in terms of found domains or enrichment in GO functions of the genes encoding them.
"""
import pandas as pd
import json
from scipy.stats import fisher_exact

"""
For a CSV file in the previously described format and two lists of identifiers (rows),
conducts a Fisher's test for each column (with Bonferroni correction) to test whether
the given subset is significantly different from the other given subset in terms of
found domains or enrichment in GO functions of the genes encoding.
"""
def perform_fishers_test(csv_file, identifiersA, identifiersB):
    df = pd.read_csv(csv_file)
    results = {}

    subsetA = df[df[df.columns[0]].isin(identifiersA)]
    subsetB = df[df[df.columns[0]].isin(identifiersB)]
    hypotheses_number = len(df.columns)
    subsetA_size = len(identifiersA)
    subsetB_size = len(identifiersB)

    for column in df.columns:
        if column == df.columns[0]: continue
        # Calculate the contingency table.
        in_subsetA = subsetA[column].sum()
        in_subsetB = subsetB[column].sum()
        table = [[in_subsetA, subsetA_size - in_subsetA],
                 [in_subsetB, subsetB_size - in_subsetB]]

        # Perform Fisher's exact test.
        _, p_value = fisher_exact(table)
        results[column] = min(1, p_value * hypotheses_number)  # With Bonferroni correction.

    result_df = pd.DataFrame(list(results.items()), columns=["ID", "p_value"])  # For easy visualisation.
    result_df_sorted = result_df.sort_values(by='p_value')

    print("Performed Fisher's test successfully.")
    return result_df_sorted


def extract_unique_values(json_filename):
    with open(json_filename, 'r') as file:
        mapping = json.load(file)
        groupA, groupB = [], []

        for queryID, list_of_seqIDs in mapping.items():
            if queryID.startswith("groupA"):
                groupA.extend(list_of_seqIDs)
            elif queryID.startswith("groupB"):
                groupB.extend(list_of_seqIDs)

    return list(set(groupA)), set(list(groupB))

if __name__ == "__main__":
    # Find the list of identifiers for group A and group B,
    # using the query_to_hits.json file.
    groupA_seqIDs, groupB_seqIDs = extract_unique_values("../../data/associations.json")

    # Check whether groups of proteins labeled as A and B in the fragments file differ significantly
    # in terms of found domains...
    # domains_result = perform_fishers_test("domains.csv", groupA_seqIDs, groupB_seqIDs)
    # print(domains_result)

    # ...and in terms of enrichment in GO functions of the genes encoding them.
    go_results = perform_fishers_test("go_terms.csv", groupA_seqIDs, groupB_seqIDs)
    print(go_results)
