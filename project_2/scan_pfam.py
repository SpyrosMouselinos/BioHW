"""
2. (6 pts) Write a program scan_pfam.py that, based on the protein file,
makes a query or queries to the hmmscan server (online server here, API description) and
downloads the result files in JSON format. Using the files downloaded from the hmmscan server,
the program should create a CSV file as output, where the rows contain the successive protein identifiers
  from the FASTA file, and the columns contain the successive PFAM protein domain identifiers.
At the intersection of a row and column, put 0 if the given domain was not found in the given protein,
   and 1 otherwise.
Unfortunately, modules for parsing HMMER output in Biopython often cause problems,
so we recommend manually reading the CSV files.
"""

import requests
import time
import pandas as pd
from Bio import SeqIO

DATA_FILE_LOCATION = "../data/"
input_fasta = DATA_FILE_LOCATION + "output_translated_proteins.fa"
output_csv = DATA_FILE_LOCATION + "pfam_domains.csv"

HMMSCAN_URL = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
HMMSCAN_RESULTS_URL = "https://www.ebi.ac.uk/Tools/hmmer/results"
HMMDB_NAME = "pfam"


def submit_hmmscan(seq):
    """
    Submit a sequence to the hmmscan service and return the job ID.
    """
    headers = {
        'Expect': '',
        'Accept': 'text/xml'
    }

    response = requests.post(url=HMMSCAN_URL,
                             headers=headers,
                             data={'hmmdb': HMMDB_NAME, 'seq': seq})
    response.raise_for_status()
    time.sleep(2)
    return response.url.split('/')[-2]


def get_results(job_id):
    """
    Retrieve the results of a completed job.
    """
    result_url = f"{HMMSCAN_RESULTS_URL}/{job_id}?output=json"
    response = requests.get(result_url)
    response.raise_for_status()
    return response.json()


def parse_results(json_result):
    """
    Parse the JSON result to extract PFAM domain hits.
    """
    domains = []
    for hit in json_result['results']['hits']:
        domains.append(hit['acc'])
    return domains


def main():
    protein_records = list(SeqIO.parse(input_fasta, "fasta"))
    protein_domains = {}

    for record in protein_records:
        seq = str(record.seq)
        protein_id = record.id
        job_id = submit_hmmscan(seq)
        json_result = get_results(job_id)
        domains = parse_results(json_result)
        protein_domains[protein_id] = domains

    # Extract unique domain identifiers
    all_domains = sorted(set(domain for domains in protein_domains.values() for domain in domains))

    # Create the 0/1 - valued CSV
    df = pd.DataFrame(0, index=protein_domains.keys(), columns=all_domains)
    for protein_id, domains in protein_domains.items():
        for domain in domains:
            df.at[protein_id, domain] = 1

    df.to_csv(output_csv)
    print(f"PFAM domains have been written to {output_csv}")


if __name__ == "__main__":
    main()
