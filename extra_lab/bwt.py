from Bio import SeqIO


class BWTGenomeSearch:
    def __init__(self, genome_input, promoter_input):
        self.genome_sequence = self.parse_input(genome_input)
        self.promoter_sequences = self.parse_input(promoter_input, is_genome=False)
        self.bwt, self.C, self.OCC = self.compute_BWT(self.genome_sequence)
        self.LF = self.last_to_first()

    def parse_input(self, input_data, is_genome=True):
        """
        Parses the input data, either from a file or directly as a string.
        """
        if input_data.endswith('.fa'):
            sequences = {}
            for record in SeqIO.parse(input_data, "fasta"):
                sequences[record.id] = str(record.seq)
            if is_genome:
                return next(iter(sequences.values()))  # Return the first sequence for genome
            else:
                return sequences  # Return all sequences for promoters
        else:
            if is_genome:
                return input_data
            else:
                return {"Test Promoter": input_data}  # Example format for direct input

    def compute_BWT(self, txt):
        """
        Computes the Burrows-Wheeler Transform (BWT) of the given text.
        Also computes the C and OCC matrices.
        """
        txt = txt + "$"
        n = len(txt)
        rotations = sorted([txt[i:] + txt[:i] for i in range(n)])
        bwt = ''.join([rotation[-1] for rotation in rotations])

        C = {}
        total_chars = sorted(set(txt))
        for char in total_chars:
            C[char] = sum(1 for c in txt if c < char)

        OCC = {}
        for char in total_chars:
            count = 0
            occ_list = []
            for c in bwt:
                if c == char:
                    count += 1
                occ_list.append(count)
            OCC[char] = occ_list

        return bwt, C, OCC

    def last_to_first(self):
        """
        Computes the Last-to-First (LF) mapping from the BWT.
        """
        LF = []
        for i, char in enumerate(self.bwt):
            LF.append(self.C[char] + self.OCC[char][i] - 1)
        return LF

    def find_occurrences(self, pattern):
        """
        Finds the occurrences of a pattern in the text using the BWT, C, and OCC matrices.
        """
        top = 0
        bottom = len(self.bwt) - 1
        while top <= bottom:
            if pattern:
                symbol = pattern[-1]
                pattern = pattern[:-1]
                if symbol in self.bwt[top:bottom + 1]:
                    top_index = self.bwt.find(symbol, top, bottom + 1)
                    bottom_index = self.bwt.rfind(symbol, top, bottom + 1)
                    top = self.C[symbol] + self.OCC[symbol][top_index] - 1
                    bottom = self.C[symbol] + self.OCC[symbol][bottom_index] - 1
                else:
                    return []
            else:
                return list(range(top, bottom + 1))
        return []

    def search_promoter_sequences(self):
        """
        Searches for promoter sequences in the genome using the BWT.
        """
        results = {}
        for promoter_id, sequence in self.promoter_sequences.items():
            occurrences = self.find_occurrences(sequence)
            results[promoter_id] = occurrences
        return results


# Test 1: Where it exists
genome_sequence = "AACGAGCCGCTTT"
promoter_sequences = "GCC"
bwt_search = BWTGenomeSearch(genome_sequence, promoter_sequences)
trivial_test_results = bwt_search.search_promoter_sequences()

print("Trivial Example Results:")
for promoter_id, positions in trivial_test_results.items():
    print(f"Promoter {promoter_id} found at positions: {positions}")

# Test 2: Where it does not exist
genome_sequence = "AACGAGCCGCTTT"
promoter_sequences = "TTA"
bwt_search = BWTGenomeSearch(genome_sequence, promoter_sequences)
trivial_test_results = bwt_search.search_promoter_sequences()

print("Trivial Example Results:")
for promoter_id, positions in trivial_test_results.items():
    print(f"Promoter {promoter_id} found at positions: {positions}")

# The actual excersize #
genome_file = "../data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa"
promoter_file = "../data/ecoli_proms.fa"
bwt_search = BWTGenomeSearch(genome_file, promoter_file)
promoter_occurrences = bwt_search.search_promoter_sequences()

for promoter_id, positions in promoter_occurrences.items():
    print(f"Promoter {promoter_id} found at positions: {positions}")
