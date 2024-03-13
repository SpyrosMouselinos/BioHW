from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Data


CODONTABLE = Data.CodonTable.unambiguous_dna_by_name["Standard"]
BLOSUM62 = substitution_matrices.load("BLOSUM62")
print(f"START Codons: {CODONTABLE.start_codons}")
print(f"STOP Codons: {CODONTABLE.stop_codons}")


def backtranslate_codons(sequence):
    return Seq(''.join([CODONTABLE.back_table[f] for f in sequence]))


def align_pair(s1, s2, mode):
    """
    Align function
    Assumption: Insertions and Deletions should only occur with "triplets" of nucleotides.
    This happens so that the ORF windows do not break while we perform an insertion or deletion.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = -1  # Use some arbitrary number
    aligner.extend_gap_score = -1  # Use a smaller penalty for extension
    aligner.substitution_matrix = BLOSUM62
    alignments = aligner.align(s1, s2)
    return alignments


def split_on_stop(text, delimiters=CODONTABLE.stop_codons):
    """
    Locates the stop condons in a sequence (on modulo 3 positions only) and returns the parts between them.
    It operates on the nucleotide level.
    It could  be much easier in the translated protein space, but I used it for debugging purposes as well.
    """
    delimiter_positions = []
    for i in range(0, len(text) - 2, 3):
        for delimiter in delimiters:
            if text[i:i + len(delimiter)] == delimiter:
                delimiter_positions.append(i)
                break
    if not delimiter_positions:
        return [text]
    pieces = []
    start = 0
    for pos in delimiter_positions:
        pieces.append(text[start:pos])
        start = pos + len(delimiters[0])
    if start < len(text):
        pieces.append(text[start:])
    return pieces


def get_all_orf(record):
    """
    Finds all open reading frames (ORFs) in a double-stranded DNA sequence.
    The logic is that we will pass the sequence twice:
     +1: Once in the positive strand
     -1: Once in the complementary strand
     Because of the nature of cloning in DNA's we must take both sequences.
     For each of them we operate in sliding windows - frames.
     We will start by not skipping any nucleotides frame: 0
     Then move on to skip 1 and then 2 (frames 1 & 2).
     There is no point in skipping 3 since this is equivalent to frame 0 with dropping the first triplet.
     So the function will scan 3 frames for 2 sequences.
     Returns the results both in nucleotide (results) and protein (translated_results) version.
    """
    translated_results = []
    results = []
    # Iterate over positive (+1) and complementary (-1) strands
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        # Check each of the three frames
        for frame in range(3):
            length = 3 * ((len(record) - frame) // 3)
            # Translate the sequence and split at stop codons which are encoded as *
            for triplet in split_on_stop(str(nuc[frame:frame + length])):
                if len(triplet) >= 3:
                    translated_results.append(Seq(triplet).translate())
                    results.append(Seq(triplet))

    return results, translated_results


def compare_sequences(s1, s2):
    """
    Compares two sequences for optimal alignment, assuming that insertions and deletions should only occur with "triplets"
     of nucleotides.
    To achieve this the function:
     1) Transforms each sequence to a set of ORFs.
     2) Picks the aligned ORF-pair with the highest score among every possible combination.
     3) Returns the codon overlap between the best pair.

    """
    results1, translated_results_1 = get_all_orf(s1)
    results2, translated_results_2 = get_all_orf(s2)
    best_score = -100
    best_seq_pair = None
    for i in translated_results_1:
        for j in translated_results_2:
            score = align_pair(i, j, 'global').score  # It is global on local subsequences
            if score >= best_score:
                best_seq_pair = (i, j)
                best_score = score

    # The optimal local alignment will be the longest common codon subsequence a.k.a the shortest of the (i,j) pair.
    optimal_alignment = best_seq_pair[0] if min([len(best_seq_pair[0]), len(best_seq_pair[1])]) == len(
        best_seq_pair[0]) else best_seq_pair[1]

    backtranslated_optimal_alignment = backtranslate_codons(optimal_alignment)
    backtranslated_optimal_alignment_t = backtranslated_optimal_alignment.reverse_complement()
    return backtranslated_optimal_alignment, backtranslated_optimal_alignment_t


# Test the provided sequences #
test_sequence_1 = "ATTGGCGAGCACA"
test_sequence_2 = "GATAGATAGCTAGCTAGTA"
answer_1 = "TTGGCGAGC"
answer_2 = "CTAGCTAGT"
s1 = SeqRecord(Seq(test_sequence_1), id="Test Seq 1")
s2 = SeqRecord(Seq(test_sequence_2), id="Test Seq 2")

print(compare_sequences(s1, s2))
