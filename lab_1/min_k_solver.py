import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import copy

### Read the Test Fasta File
sequences = []
ids = []
lengths = []
for seq_record in SeqIO.parse("../data/yeast.fa", "fasta"):
    id = seq_record.id
    seq = str(seq_record.seq)
    record_len = len(seq_record)
    sequences.append(seq)
    ids.append(id)
    lengths.append(record_len)

max_possible_k_value = min(lengths)
print(max_possible_k_value)


### Are there any two sequences that are exactly the same?
def check_hash_match(strings):
    seen = set()
    for string in strings:
        string_hash = hash(string)
        if string_hash in seen:
            return True
        seen.add(string_hash)
    return False


# Check for hash match
#hash_match_found = check_hash_match(sequences)
# print(hash_match_found)


def kmers(s, k):
    """
    Generates a set of k-mers
    :param s: The sequence
    :param k: The value of k
    :return: A set of k-mers
    """
    kmers = set()
    for i in range(len(s) - k + 1):
        kmer = s[i:i + k]
        if kmer in kmers:
            pass
        else:
            kmers.add(kmer)
    return kmers


sequence_hashmap = {}
kmer_hashmap = {}

current_k = 464 # Correct Answer found after binary search over 2-700
# For each new sequence
for sid, sequence in enumerate(sequences):
    sequence_hashmap[sid] = set()
    # For k-mer in sequence
    for kmer in kmers(sequence, current_k):
        # Have we seen it so far?
        if kmer in kmer_hashmap:
            # If yes, go over the sequences we have seen it over and delete it from their sets.
            for seen_sequence in kmer_hashmap[kmer]:
                if seen_sequence != sid:
                    # Go to sequence hashmap and delete it from the sequence
                    sequence_hashmap[seen_sequence].remove(kmer)
                    # Check if you left any sequence empty
                    if len(sequence_hashmap[seen_sequence]) == 0:
                        print("Not possible")
                        sys.exit(1)
                else:
                    pass  # If we have seen it in the same sequence, do not do anything
            # Remove it from the kmer_hashmap now so not to do double work
            kmer_hashmap[kmer] = []
        else:
            # If we have not seen it so far, update the two hashmaps.
            kmer_hashmap[kmer] = [sid]
            sequence_hashmap[sid].update([kmer])

    # At the end, check if anything was added to the sequence_hashmap
    if len(sequence_hashmap[sid]) == 0:
        print("Not possible")
        sys.exit(1)

print(f"Possible at {current_k}!")