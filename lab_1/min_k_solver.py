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
hash_match_found = check_hash_match(sequences)
if hash_match_found:
    print("There are at least 2 same sequences in the set")
else:
    print("All sequences are unique!")


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


# Solution Start #


sequence_hashmap = {}
kmer_hashmap = {}




def solve(sequences, current_k):
    """
    This function will act as the main logic to solve the problem. I make the following assumptions to make things
    faster.
    # Assumption 1: Since we are only interested in k, we do not need to find the complementary probes but just the
        unique subsequences per sequence. Thus instead of performing .reverse_complement() to each sequence we can omit
         it.
    # Assumption 2: We can bound k in the range [1, min(len(sequence))], since any k larger than the shortest sequence
        will not solve the problem. Furthermore, since we operate on unique sequences (we checked above), we can create
         at least log_4(len(sequences)) without overlap ~ log_4(4437) = 6. So the final range is [6, 702]
    # Assumption 3: While we can not make the solution faster, we can quickly find out when a current execution will not
        lead to a correct result. It's very simple: If for any sequence, we have so far seen all k-mers that exist in
        it,there is no point in executing more steps since there is at least 1 sequence that will not have a unique
        probe. This is why we maintain a hashmap per sequence where we keep track of all the so far seen k-mers, if
        we find that all the k-mers of any sequence have been seen, we terminate.
    # Assumption 4: Apart from the sequence hashmap we can maintain another hashmap on the k-mers, to check if we have
        already seen that k-mer. If yes, we go over all sequences where it exists and mark it as seen.
    """

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
                            return False
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
            return False

    print(f"Possible at {current_k}!")
    return True


search_range = [6, 702]


def binary_search(low, high):
    print(f"Checking k={high}\n")
    if not solve(sequences, high):
        return "No solution found within the given range."
    while low < high:
        mid = (low + high) // 2
        print(f"Checking k={mid}\n")
        if solve(sequences, mid):
            high = mid
        else:
            low = mid + 1
    return low


result = binary_search(6, 702)
print(result)

#k = 464 Correct Answer