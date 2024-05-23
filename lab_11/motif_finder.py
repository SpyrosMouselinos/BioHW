import numpy as np
import random
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.motifs import write

class MotifFinder:
    def __init__(self, sequences, k):
        """
        Initialize the MotifFinder with sequences and motif length.

        Args:
        sequences (list of SeqRecord): List of DNA sequences.
        k (int): Length of the motif.
        """
        self.sequences = [str(seq.seq) for seq in sequences]
        self.k = k
        self.nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.bg = self.calculate_background_probabilities(k)

    def add_pseudocounts(self, pspm, pseudocount=0.035):
        """
        Add pseudocounts to the PSPM to avoid zero probabilities.
        Returns PSPM with added pseudocounts.
        Important note:
        I finetuned the pseudocount value so that I get information_content scores similar to those in the examples in the lecture slides.
        Alternatively i could have used the SQRT(N) * BG where N is the depth of the iteration
        """
        return pspm + pseudocount

    def normalize_pspm(self, pspm):
        """
        Normalize the PSPM so that each column sums to 1.
        We assume A C G T as column items and position as row
        Normalized PSPM.
        """
        return pspm / pspm.sum(axis=0, keepdims=True)

    def calculate_pwm(self, pspm):
        """
        Calculate the PWM from the PSPM and background probabilities.
        Returns the Position Weight Matrix.
        """
        return np.log2(pspm / self.bg)

    def calculate_information_content(self, pspm):
        """
        Calculate the information content for each position in the PSPM.
        Returns the Total information content.
        """
        ic = 0  # For debugging, we dont need an array.
        for i in range(pspm.shape[0]):
            for j in range(pspm.shape[1]):
                ic += pspm[i, j] * np.log2(pspm[i, j] / self.bg[i, j])
        return ic

    def score_sequence(self, sequence, pwm):
        """
        Score a sequence against the PWM.
        Returns the Score of the sequence based on the PWM.
        """
        score = 0.0
        for i, nucleotide in enumerate(sequence):
            score += pwm[self.nucleotide_to_index[nucleotide], i]
        return score

    def build_pspm(self, sequences):
        """
        Build the PSPM from the given sequences.
        Returns the Position Specific Probability Matrix.
        """
        pspm = np.zeros((len(self.nucleotide_to_index), self.k))
        for seq in sequences:
            for i, nucleotide in enumerate(seq):
                pspm[self.nucleotide_to_index[nucleotide], i] += 1
        return pspm

    def calculate_background_probabilities(self, k):
        """
        Calculate background probabilities from all given sequences.
        """
        bg_probabilities = np.ones(shape=(4, k)) * 0.25
        return bg_probabilities

    def greedy_consensus_algorithm(self):
        """
        Perform the greedy consensus algorithm to find motifs.
        Returns the List of top motifs with their information content and sequences.
        """
        best_motifs = []
        first_sequence = self.sequences[0]

        # Iterate over each starting position in the first sequence
        for i in range(len(first_sequence) - self.k + 1):
            initial_motif = first_sequence[i:i + self.k]
            pspm = self.build_pspm([initial_motif])
            pspm = self.add_pseudocounts(pspm)
            pspm = self.normalize_pspm(pspm)
            # print(self.calculate_information_content(pspm)) # Debug
            selected_sequences = [initial_motif]

            # Iterate over each subsequent sequence
            for seq in self.sequences[1:]:
                best_subseq = None
                best_score = -np.inf

                # Find the best fitting subsequence
                for j in range(len(seq) - self.k + 1):
                    subseq = seq[j:j + self.k]
                    pwm = self.calculate_pwm(pspm)
                    score = self.score_sequence(subseq, pwm)

                    if score > best_score:
                        best_score = score
                        best_subseq = subseq

                selected_sequences.append(best_subseq)
                pspm = self.build_pspm(selected_sequences)
                pspm = self.add_pseudocounts(pspm)
                pspm = self.normalize_pspm(pspm)

            ic = self.calculate_information_content(pspm)
            best_motifs.append((ic, selected_sequences))

        # Sort motifs by information content and return the top 5
        best_motifs.sort(reverse=True, key=lambda x: x[0])
        return best_motifs[:5]


    def save_motifs_to_pfm(self, sequences, filename):
        """
        Save the motifs to a file in .pfm format suitable for the JASPAR database.
        sequences: List of motifs sequences.
        filename: The filename to save the motifs.
        """
        for id, motif_group in enumerate(sequences):
            motif_list = motifs.create(motif_group[0]).format('pfm')
            with open(f'output_motif_{id}.pfm', 'w') as fout:
                fout.write(motif_list)
        return


def run_multiple_motif_searches(sequences, k, num_runs):
    """
    Run the greedy consensus algorithm multiple times with shuffled sequences.

    Args:
    sequences (list of SeqRecord): List of sequences.
    k (int): Length of the motif.
    num_runs (int): Number of times to run the algorithm with shuffled sequences.

    Returns all results and the top 5 most common motifs.
    """
    all_results = []

    for _ in range(num_runs):
        shuffled_sequences = sequences.copy()
        random.shuffle(shuffled_sequences)
        motif_finder = MotifFinder(shuffled_sequences, k)
        motifs = motif_finder.greedy_consensus_algorithm()
        all_results.extend(motifs)

    # Count the occurrences of each motif sequence
    motif_counter = Counter(tuple(motif[1]) for motif in all_results)
    most_common_motifs = motif_counter.most_common(5)

    return all_results, most_common_motifs


def test_sanity():
    """
    Run a sanity test with predefined sequences from lecture slides.

    Sequences:
    Sequence 1: ACTGA
    Sequence 2: TAGCG
    Sequence 3: CTTGC

    Motif length: 4
    Runs the algorithm once with equal probability of ACTG symbols.
    """
    sequences = [
        SeqRecord(Seq("ACTGA"), id="seq1"),
        SeqRecord(Seq("TAGCG"), id="seq2"),
        SeqRecord(Seq("CTTGC"), id="seq3")
    ]
    k = 4

    motif_finder = MotifFinder(sequences, k)
    motifs = motif_finder.greedy_consensus_algorithm()

    # Print all results
    print("Sanity Test Results:")
    for ic, motif in motifs:
        print(f"Information Content: {ic}")
        print(f"Motif Sequences: {motif}")


def main():
    # Read the sequences from the file
    sequences = list(SeqIO.parse("../data/motif_data.fa", "fasta"))
    k = 7
    num_runs = 1

    # Run the motif search multiple times
    all_results, most_common_motifs = run_multiple_motif_searches(sequences, k, num_runs)

    # Save the results to a .pfm file
    motif_finder = MotifFinder(sequences, k)
    motif_finder.save_motifs_to_pfm(most_common_motifs, "pfm")

    # Print all results
    print("All Results:")
    for ic, motif in all_results:
        print(f"Information Content: {ic}")
        print(f"Motif Sequences: {motif}")



if __name__ == "__main__":
    main()

    # By running my 5 top motifs against the JASPAR Database online I get hits so i assume everything runs correctly #

