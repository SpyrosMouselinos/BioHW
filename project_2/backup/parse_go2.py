
"""
Assigns Gene Ontology functions to genes found based on fragments (based on the e_coli.gaf file)
and returns the result in a CSV file in the format : gene identifiers in rows, GO identifiers in columns,
1 at the intersection if the function is assigned to the gene, 0 otherwise.
"""
import csv
from Bio import SeqIO

"""
Extracting gene ids from the fasta file.
"""
def parse_fasta(fasta_file):
    gene_ids = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_ids.add(record.id)
    return gene_ids

"""
Based on the .gaf file and gene identifiers, returns mapping of the gene ids to
the genes' GO function, as well as set of all GO ids.
"""
def parse_gaf(gaf_file, gene_ids):
    gene_to_go_mapping = {}
    all_go_terms = set()

    with open(gaf_file, 'r') as file:
        _ = file.readline()
        for line in file:
            columns = line.strip().split('\t')
            gene_id = columns[2]  # The gene ID
            go_id = columns[4]  # GO ID

            if gene_id in gene_ids:
                if gene_id not in gene_to_go_mapping.keys():
                    gene_to_go_mapping[gene_id] = set()
                gene_to_go_mapping[gene_id].add(go_id)
                all_go_terms.add(go_id)

    # The genes for which no function was assigned.
    for gene_id in gene_ids:
        if gene_id not in gene_to_go_mapping.keys():
            gene_to_go_mapping[gene_id] = set()

    return gene_to_go_mapping, all_go_terms


"""
Creates a CSV file in the format : gene identifiers in rows, GO identifiers in columns,
1 at the intersection if the function is assigned to the gene, 0 otherwise.
"""
def write_csv(gene_to_go_mapping, all_go_terms, output_filename):
    with open(output_filename, 'w', newline='') as csvfile:
        column_names = ['Seq ID'] + sorted(all_go_terms)
        writer = csv.DictWriter(csvfile, fieldnames=column_names)
        writer.writeheader()

        for gene_id, gene_go_terms in gene_to_go_mapping.items():
            row = {'Seq ID' : gene_id}
            for go_term in all_go_terms:
                row[go_term] = '1' if go_term in gene_go_terms else '0'
            writer.writerow(row)
    print("Successfully created the .csv file.")

if __name__ == "__main__":

    # Extracting names of the sequences from database that were hit.
    gene_ids = parse_fasta("hit_proteins.fa")

    # Creating a mapping from those genes into their function and creating the output csv file.
    gene_to_go_mapping, all_go_terms = parse_gaf("../../data/e_coli.gaf", gene_ids)
    write_csv(gene_to_go_mapping, all_go_terms, "go_terms.csv")
