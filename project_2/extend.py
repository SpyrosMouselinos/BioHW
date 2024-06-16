"""
1. (6 pts) Write a program extend.py that,
for a given list of protein fragments in FASTA format (protein_fragments.fa),
finds the closest protein in a local database created from a given larger FASTA file (genes_e_coli.fa)
and returns a new FASTA file with the resulting gene sequences already translated into proteins.
Hint: how to run BLAST locally.
Please note that we have fragments of protein amino acid sequences and DNA sequences of genes.
"""
import json
import subprocess
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

DATA_FILE_LOCATION = "../data/"
protein_fragments_file = DATA_FILE_LOCATION + "protein_fragments.fa"
genes_e_coli_file = DATA_FILE_LOCATION + "genes_e_coli_new.fa"
genes_e_coli_db = "DB"
output_file = DATA_FILE_LOCATION + "output_translated_proteins.fa"
associations_filename = DATA_FILE_LOCATION + "associations.json"

# Inspect the files #
from Bio import SeqIO


def inspect_fasta_file(fasta_file, num_sequences=3):
    """
    Inspect the content of a FASTA file.
    """

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"File: {fasta_file}")
    print(f"Total sequences: {len(sequences)}\n")

    for i, seq_record in enumerate(sequences[:num_sequences], 1):
        print(f"Sequence {i}:")
        print(f"ID: {seq_record.id}")
        print(f"Description: {seq_record.description}")
        print(f"Sequence length: {len(seq_record.seq)}")
        print(f"Sequence: {seq_record.seq[:50]}...")
        print("-" * 40)


### MAIN EXCERSIZE STARTS ###
# Create Nucleotide BLAST database from the e_coli
subprocess.run(f"makeblastdb -in {genes_e_coli_file} -out {genes_e_coli_db} -dbtype nucl -title ECOLI_NUCL_DB",
               shell=True)

# Read protein fragments
protein_fragments = list(SeqIO.parse(protein_fragments_file, "fasta"))
genes_e_coli = SeqIO.to_dict(SeqIO.parse(genes_e_coli_file, "fasta"))

# Prepare the output list for translated protein sequences
translated_proteins = []

### Create an associations file between groups and ids ###
associations = {}

# Loop over each protein fragment and run BLAST
for fragment in protein_fragments:
    associations.update({fragment.id: []})
    fragment_file = "temp_fragment.fa"
    with open(fragment_file, "w") as temp_f:
        SeqIO.write(fragment, temp_f, "fasta")

    # Run TBLASTN to find the closest match in the nucleotide database
    # We run TBLASTN because we run Proteins against Nucleotides
    blast_output = "temp_blast.xml"
    blast_cmd = f"tblastn -query {fragment_file} -db C:\\Users\\mouse\\PycharmProjects\\BioHW\\project_2\\{genes_e_coli_db} -out {blast_output} -evalue 0.001 -outfmt 5"
    subprocess.run(blast_cmd, shell=True)

    # Parse the BLAST output
    blast_records = list(SearchIO.parse(blast_output, "blast-xml"))

    # # Get the best hit for the current fragment
    if blast_records and blast_records[0].hits:
        for hit in range(len(blast_records[0].hits)):
            protein_id = blast_records[0].hits[hit].id
            # Retrieve the corresponding Translated GenSeq into Protein
            translated_gene_sequence_into_protein = genes_e_coli[protein_id].seq.translate()

            # Rename some stuff for clarity
            translated_protein = SeqRecord(str(translated_gene_sequence_into_protein),
                                           description="",
                                           id=protein_id)

            # Add the translated protein to the output list
            translated_proteins.append(translated_protein)
            associations[fragment.id].append(protein_id)

with open(associations_filename, "w") as out:
    json.dump(associations, out)

# Write the translated protein sequences to the output FASTA file
with open(output_file, "w") as output_f:
    SeqIO.write(translated_proteins, output_f, "fasta")

print(f"Translated hits protein sequences have been written to {output_file}")
