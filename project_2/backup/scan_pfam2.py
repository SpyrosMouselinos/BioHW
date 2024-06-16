import requests
from Bio import SeqIO
import json
import csv
import time
"""
Takes sequence as plain text and performs hmmscan via website
(searches a profile HMM against a sequence database).
"""
def perform_hmmscan(sequence, filename):
    url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
    headers = {'Expect': '',
               'Accept': 'application/json'
    }
    data = {
        'hmmdb': (None, 'pfam'),
        'seq': (None, sequence)
    }
    # Send the request
    response = requests.post(url, headers=headers, files=data)

    # Write the result to a json file.
    if response.status_code == 200:
        with open(filename, 'w') as file:
            file.write(response.text)
        print("hmmscan performed successfully; result written into file:", filename)
    else:
        print("Failed to perform hmmscan. Status code :", response.status_code)
        print("Response :", response.text)

"""
Takse list of json_files' names, parses the results into a dictionary
of sequence names as keys and set of hit domains as values for each sequence's name.
"""
def parse_results(json_files, protein_ids):
    found_domains = {}
    all_domains = set()

    for json_file, protein_id in zip(json_files, protein_ids):
        print("Parsing file : ", json_file)
        with open(json_file, 'r', encoding='utf-8') as file:
            json_data = json.load(file)
            for result in json_data['results']['hits']:
                hit_domains = set(domain['alihmmname'] for domain in result['domains'])

                found_domains[protein_id] = hit_domains
                all_domains.update(hit_domains)

    return found_domains, all_domains

"""
Creates a csv file, from a list of json files, where the rows contain the successive
protein identifiers and the columns contain the successive PFAM protein domain identifiers.
At the intersection of a row and column, puts 0 if the given domain was not found in
the given protein, and 1 otherwise. 
"""
def create_csv(found_domains, all_domains, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        column_names = ['Seq ID'] + sorted(all_domains)
        writer = csv.DictWriter(csvfile, fieldnames=column_names)
        writer.writeheader()

        for protein, protein_domains in found_domains.items():
            row = {'Seq ID': protein}
            for domain in all_domains:
                row[domain] = '1' if domain in protein_domains else '0'
            writer.writerow(row)
    print("Successfully created the .csv file.")

if __name__ == "__main__":
    proteins_list = list(SeqIO.parse("../data/output_translated_proteins.fa", "fasta"))

    # Getting the hmmscans for the protein sequences.
    json_files = []
    protein_ids = []
    for record in proteins_list:
        filename = f"json_result_{record.id}.json"
        #perform_hmmscan(str(record.seq), filename)
        #time.sleep(1)
        json_files.append(filename)
        protein_ids.append(record.id)

    # Parsing obtained results and creating the csv file.
    found_domains, all_domains = parse_results(json_files, protein_ids)
    create_csv(found_domains, all_domains, "domains.csv")
