"""
3. (4 pts) Write a program parse_go.py that assigns Gene Ontology functions to genes found based on fragments
(based on the e_coli.gaf file) and returns the result in a CSV file in the format as in point 2
(gene identifiers in rows, GO identifiers in columns, 1 at the intersection if the function is assigned to the gene,
0 otherwise).
"""

import pandas as pd
from Bio import SeqIO
from Bio.UniProt.GOA import gafiterator

DATA_FILE_LOCATION = "../data/"
translated_proteins_file = DATA_FILE_LOCATION + "output_translated_proteins.fa"
gaf_file = DATA_FILE_LOCATION + "e_coli.gaf"
output_csv = DATA_FILE_LOCATION + "go_annotations.csv"


def fix_gaf_file(file):
    flag = False
    with open(file, 'r') as fin:
        inline = fin.readline()
        if inline.strip() != "!gaf-version: 2.0":
            flag = True

    if flag:
        with open(file, 'r') as fin:
            store_lines = fin.readlines()

        # Add the custom string as the first line
        store_lines.insert(0, '!gaf-version: 2.0\n')

        with open(file, 'w') as fout:
            fout.writelines(store_lines)
    return file


def parse_gaf(file):
    """
    Parse the GAF file to extract gene and GO term associations.
    Note: I met a problem reading the file because there was no GAF version in the first line so I add it manually.
    From what I saw online our file fits the 2.0 version not the 1.0 version
    """
    gene_go_mapping = {}

    with open(file, 'r') as handle:
        for rec in gafiterator(handle):
            gene_id = rec['DB_Object_Symbol']
            go_id = rec['GO_ID']
            if gene_id not in gene_go_mapping:
                gene_go_mapping[gene_id] = set()
            gene_go_mapping[gene_id].add(go_id)

    return gene_go_mapping


def get_gene_ids_from_fasta(fasta_file):
    """
    Extract gene identifiers from the translated protein FASTA file.
    """
    gene_ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        description_parts = record.description.split(" ")
        gene_id = description_parts[-1]  # Assuming the gene ID is the last part (That's what I did in extend.py)
        gene_ids.append(gene_id)
    return gene_ids


def create_go_matrix(genes, gene_go_mapping):
    """
    Create a DataFrame where rows are genes and columns are GO terms.
    """
    all_genes = sorted(set(genes))
    all_go_terms = sorted({go for gos in gene_go_mapping.values() for go in gos})

    df = pd.DataFrame(0, index=all_genes, columns=all_go_terms)

    for gene in genes:
        if gene in gene_go_mapping:
            for go in gene_go_mapping[gene]:
                df.at[gene, go] = 1

    return df


def main():
    genes = get_gene_ids_from_fasta(translated_proteins_file)
    gene_go_mapping = parse_gaf(fix_gaf_file(gaf_file))
    go_matrix = create_go_matrix(genes, gene_go_mapping)
    go_matrix.to_csv(output_csv)
    print(f"GO annotations have been written to {output_csv}")


if __name__ == "__main__":
    main()
