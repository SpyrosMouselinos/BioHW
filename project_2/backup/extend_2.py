import subprocess
import json
from Bio import SeqIO

"""
Compare query protein sequences against a nucleotide database,
using local BLAST.
"""


def run_tblastn_locally(protein_query, nucleotide_db, output_file):
    tblastn_command = [
        "tblastn",
        "-query", protein_query,
        "-db", nucleotide_db,
        "-evalue", "0.001",  # Treshold for the E-values
        "-outfmt", "6",  # Format of the output.
        "-out", output_file
    ]

    result = subprocess.run(tblastn_command, capture_output=True, text=True)
    if result.returncode != 0:
        print("BLAST error:", result.stderr)
    else:
        print("BLAST completed successfully.")




# def create_database(nucleotide_seq, nucleotide_db):
#     makeblastdb_command = [
#         "makeblastdb",
#         "-in", nucleotide_seq,
#         "-dbtype", "nucl",
#         "-out", nucleotide_db
#     ]
#     subprocess.run(makeblastdb_command, check=True)


"""
Based on .txt BLAST's output file, extracts all the hit protein sequences,
and writes it to a fasta file. Additionaly, maps the identifiers of the query sequences
to the lists of their corresponding hits ids, and writes this mapping to a json file.
"""


def extract_proteins(blast_output_file, nucleotide_seq, output_filename, json_filename):
    hit_positions = []
    query_to_hits = {}

    with open(blast_output_file) as blast_file:
        for line in blast_file:
            cols = line.split()

            query_seq_id = cols[0]
            hit_protein_id = cols[1]
            hit_positions.append(hit_protein_id)
            if query_seq_id not in query_to_hits.keys():
                query_to_hits[query_seq_id] = []
            query_to_hits[query_seq_id].append(hit_protein_id)

    records = SeqIO.to_dict(SeqIO.parse(nucleotide_seq, "fasta"))
    hit_proteins = []

    for protein_id in set(hit_positions):
        records[protein_id].seq = records[protein_id].seq.translate()
        hit_proteins.append(records[protein_id])

    with open(output_filename, "w") as output_file:
        SeqIO.write(hit_proteins, output_file, "fasta")
    print("Extracting the hit proteins completed successfully. The output saved to ", output_filename)

    with open(json_filename, "w") as json_file:
        json.dump(query_to_hits, json_file)
    print("Saving the additional query's id-to-hit mapping completed successfully. The ouput saved to ", json_filename)


if __name__ == "__main__":
    #Running the BLAST on given protein fragments.
    run_tblastn_locally(
        "../data/protein_fragments.fa",
        "DB",
        "blast_results.txt"
    )

    # Extracting the hit proteins into a fasta file.
    extract_proteins(
        "blast_results.txt",
        "../data/genes_e_coli_new.fa",
        "hit_proteins.fa",
        "query_to_hits.json"
    )
