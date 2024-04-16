import numpy as np
from Bio import Align
from Bio import Phylo
from Bio.Phylo.TreeConstruction import  DistanceTreeConstructor, DistanceMatrix
from Bio.Seq import Seq
from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")


class Profile:
    """
    A simple class representing a profile which is one or more sequences.
    """

    def __init__(self, sequences):
        # Small Hack for calculating the DP easier: Append a space '-' at the beginning of each sequence
        self.sequences = ['-' + f for f in sequences]
        self.alignment_length = len(self.sequences[0])
        for s in self.sequences:
            assert len(s) == self.alignment_length

    def shape(self):
        return len(self.sequences), self.alignment_length

    def get_sequences(self):
        return [Seq(f) for f in self.sequences]

    def add_sequence(self, new_sequence):
        assert len(new_sequence) == self.alignment_length
        self.sequences.append(new_sequence)

    def get_column(self, index):
        assert index >= 0
        assert index <= self.alignment_length
        return [f[index] for f in self.sequences]

    def pprint(self):
        # Pretty print the alignment #
        print("MSA Result:\n")
        for i in range(len(self.sequences)):
            print(f"Sequence {i}: {self.sequences[i][1:]}\n")


def generate_guide_upgma_tree(sequences, sub_mat, gap_score=-4):
    """
    Generates the guide tree from calculating the distance matrix over all sequence pairs.
    The method UPGMA is used to generate the tree.
    :param sequences: The sequences to be aligned
    :param sub_mat: The substitution matrix to be used.
    :return:
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load(sub_mat)
    aligner.gap_score = gap_score

    # Calculation of a Distance Matrix
    # Here, I am going use a very simple method.
    # I know that in the class slides Clustal uses the Kimura distance but
    # Since this is a simple demonstration, for each pairwise alignment, i penalize by 1 for every mismatch and 0.5 per gap
    # For example AAA and AAA return 0
    # AA- and AAB return 0.5
    # while AA- and BAB return 1.5
    def simple_distance_of_sequences(seq_1, seq_2):
        mismatch = 0
        for a, b in zip(seq_1, seq_2):
            if a == '-' and b != '-':
                mismatch += 0.5
            elif a != '-' and b == '-':
                mismatch += 0.5
            else:
                if a != b:
                    mismatch += 1
        return mismatch / len(seq_1)

    distance_matrix_lower_triangular_format = []

    for i, seq_1 in enumerate(sequences):
        distance_matrix_row = []
        for j, seq_2 in enumerate(sequences[:(i + 1)]):
            alignment = aligner.align(Seq(seq_1), Seq(seq_2))[0]
            aligned_sequence_1 = alignment[0]
            aligned_sequence_2 = alignment[1]
            distance_matrix_row.append(simple_distance_of_sequences(aligned_sequence_1, aligned_sequence_2))
        distance_matrix_lower_triangular_format.append(distance_matrix_row)

    bio_distance_matrix = DistanceMatrix(sequences, distance_matrix_lower_triangular_format)

    # Construct the UPGMA tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(bio_distance_matrix)
    # Print the guide tree
    Phylo.draw_ascii(tree)
    return tree


def calculate_sp_score(profile_1_items, profile_2_items, sub_mat=None, mismatch_penalty=-3, gap_penalty=-2):
    """
    Function that calculates the SP score between two sequences or profiles to be compared.
    This is a simple column operation that compares each item of profile_1 with each item of profile_2
    Use the provided score matrix to update the substitution scores.
    """
    total_score = 0
    # Calculate SP score by summing up all pairs between two columns
    for element1 in profile_1_items:
        for element2 in profile_2_items:
            if sub_mat is not None:
                if (element1, element2) in sub_mat:
                    total_score += sub_mat[(element1, element2)]
                elif element1 == '-' and element2 == '-':
                    total_score += 0  # Gap to Gap penalty
                else:
                    total_score += gap_penalty  # Item to Gap penalty
            else:
                if element1 == element2:
                    total_score += 0
                elif element1 == '-' and element2 != element1:
                    total_score += gap_penalty  # Item to Gap penalty
                elif element2 == '-' and element2 != element1:
                    total_score += gap_penalty  # Item to Gap penalty
                else:
                    total_score += mismatch_penalty
    return total_score


def profile_align(profile1, profile2, sub_mat=None, mismatch_penalty=-3, gap_penalty=-2):
    """
    Dynamic programming construction matrix comparing two profiles containing a variable number of sequences.
    """
    # Initialize the DP matrix
    M = profile2.shape()[1]
    N = profile1.shape()[1]
    dp = np.zeros(shape=(M, N))
    for i in range(N):
        dp[0, i] = i * gap_penalty * max(profile1.shape()[0], profile2.shape()[0])
    for j in range(M):
        dp[j, 0] = j * gap_penalty * max(profile1.shape()[0], profile2.shape()[0])

    # Fill the DP matrix
    for i in range(1, M):
        for j in range(1, N):
            s_diag = calculate_sp_score(profile1.get_column(j), profile2.get_column(i),
                                        sub_mat=sub_mat,
                                        mismatch_penalty=mismatch_penalty, gap_penalty=gap_penalty)
            score_diagonal = dp[i - 1][j - 1] + s_diag
            score_up = dp[i - 1][j] + gap_penalty * max(profile1.shape()[0], profile2.shape()[0])
            score_left = dp[i][j - 1] + gap_penalty * max(profile1.shape()[0], profile2.shape()[0])
            dp[i][j] = max(score_diagonal, score_up, score_left)

    return dp


def backtrack_alignment(dp, profile1, profile2, sub_mat=None, mismatch_penalty=-3, gap_penalty=-2):
    """
    Backtracking algorithm on the results of the profile_align dp matrix.
    Returns the joint aligned profile as a result of aligning profile1 and profile2
    """
    i, j = profile2.shape()[1] - 1, profile1.shape()[1] - 1
    alignment1, alignment2 = [], []

    # Traverse the DP matrix from bottom-right to top-left
    while i > 0 and j > 0:
        if dp[i][j] == dp[i - 1][j - 1] + calculate_sp_score(profile1.get_column(j), profile2.get_column(i),
                                                             sub_mat=sub_mat,
                                                             mismatch_penalty=mismatch_penalty,
                                                             gap_penalty=gap_penalty):
            alignment1.insert(0, profile1.get_column(j))
            alignment2.insert(0, profile2.get_column(i))
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty * max(profile1.shape()[0], profile2.shape()[0]):
            alignment1.insert(0, profile1.get_column(j - 1))
            alignment2.insert(0, ['-'])
            i -= 1
        else:
            alignment1.insert(0, ['-'])
            alignment2.insert(0, profile2.get_column(i - 1))
            j -= 1

    # Handle the remaining columns
    while i > 0:
        alignment1.insert(0, profile1.get_column(j - 1))
        alignment2.insert(0, ['-'])
        i -= 1

    while j > 0:
        alignment1.insert(0, ['-'])
        alignment2.insert(0, profile2.get_column(i - 1))
        j -= 1

    # Reconstruct profiles from alignments #
    final_alignment = []
    for a, b in zip(alignment1, alignment2):
        final_alignment.append(a + b)

    new_aligned_profile = Profile([''.join(group) for group in zip(*final_alignment)])
    return new_aligned_profile


def progressive_alignment(tree,
                          sub_mat=None,
                          mismatch_penalty=-3,
                          gap_penalty=-2):
    """
    A DFS algorithm (recursive) that visits the guide tree and performs profile alignments progressively.
    It assumes that the guide tree is bifurcated (UPGMA creates such trees for example)
    """
    # Prepare a dictionary to store the profiles/alignments by clade
    alignments = {leaf.name: Profile([leaf.name]) for leaf in tree.get_terminals()}

    # Recursive function to align profiles
    def align_profiles(node):
        # Base case, node is a leaf --> Return the leaf as a profile #
        if node.is_terminal():
            return alignments[node.name]

        # Node is not leaf --> Align left with right children and return a new profile #
        left, right = node.clades
        left_profile = align_profiles(left)
        right_profile = align_profiles(right)
        dp_matrix = profile_align(left_profile,
                                  right_profile,
                                  sub_mat=sub_mat,
                                  mismatch_penalty=mismatch_penalty,
                                  gap_penalty=gap_penalty)

        # Since profile_align returns DP matrix, use backtrack to get alignment
        aligned = backtrack_alignment(dp_matrix, left_profile,
                                      right_profile,
                                      sub_mat=sub_mat,
                                      mismatch_penalty=mismatch_penalty,
                                      gap_penalty=gap_penalty)
        alignments[node.name] = aligned
        return alignments[node.name]

    # Start alignment from the root
    final_alignment = align_profiles(tree.root)
    return final_alignment


# Very easy example that can be reproduced below #
sequences = ["ACTCAT", "AGTCAT", "ACGTCCT"]

# Step 1: Load the guide tree and the Sequences.
# For this code I also included a function to use the UPGMA to generate the guide tree, so that I can debug
# my code.

# To read a tree from a file we would have used the code below
# guide_tree = Phylo.read("guide_tree.nw", "newick")

# Or with my code
guide_tree = generate_guide_upgma_tree(sequences, sub_mat='BLOSUM62')

# Step 2:
# Use the guide tree to perform the MSA recursively
final_alignment = progressive_alignment(guide_tree, sub_mat=None, mismatch_penalty=-3, gap_penalty=-2)
final_alignment.pprint()
