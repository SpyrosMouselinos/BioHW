# Project 1
# Mouselinos Spyridon

# Imports
from Bio import SeqIO
from Bio import Data
from Bio.Seq import Seq
from Bio import Phylo
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
import subprocess
from Bio import AlignIO

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator

# Step 1:
# Read the given protein sequence from a .fa file and generate one
# of the possible DNA sequences that could encode it
# (a proper codon can be randomly selected for each amino acid).


PROTEIN_SEQUENCE_FILE = '../data/TPH2.fa'
TREE_FILE = '../data/tree'
CODONTABLE = Data.CodonTable.unambiguous_dna_by_name["Standard"]

for seq_record in SeqIO.parse(PROTEIN_SEQUENCE_FILE, "fasta"):
    id = seq_record.id
    original_sequence = str(seq_record.seq)
    record_len = len(seq_record)


def backtranslate_codons(sequence):
    """
    Selects a single CODON to backtranslate a protein
    sequence to base sequence.
    """
    return Seq(''.join([CODONTABLE.back_table[f] for f in sequence]))


back_translated_sequence = backtranslate_codons(original_sequence)
print(
    f"Original sequence:\n{original_sequence[:5]}...\nOriginal Length: {record_len}\n\n"
    f"Backtranslated sequence:\n{back_translated_sequence[:5]}...\n"
    f"Backtranslated Length: {len(back_translated_sequence)}")

# Step 2:
# Load the tree data from a file in Newick format
# (for our tree, its description in Newick format looks as it does in the "tree" file).

original_tree = Phylo.read(TREE_FILE, 'newick')
# Let's visualise it
Phylo.draw(original_tree, branch_labels=lambda c: c.branch_length)


# Step 3:
# Randomly generate sequences A-E using a discrete-time Markov model with assumptions similar to the Jukes-Cantor model
# (the probability of mutation in one step of the simulation should be 1/5000 per position).
# Choose the number of steps in the algorithm so that the expected number of mutations matches the given tree.

def markov_jukes_cantor(seq, t, mu=1.0 / 5000):
    """
    Since we assume that the probability of mutation in one step (t) is equal to 1/5000 per position
    we can set mu = 1/5000.

    :param seq: The initial sequence to mutate
    :param t: Timestep
    :param mu: Mu parameter in P-matrix of JC model
    """

    def calculate_p_matrix(t, mu=1.0 / 5000):
        """Calculate the transition probability matrix P for the Jukes-Cantor model."""
        e_factor = np.exp(-4.0 * mu * t / 3.0)
        p_matrix = np.full((4, 4), (1.0 - e_factor) / 4.0)
        np.fill_diagonal(p_matrix, 1.0 / 4.0 + 3.0 / 4.0 * e_factor)
        return p_matrix

    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    index_to_nucleotide = {index: nt for nt, index in nucleotide_to_index.items()}
    evolved_sequence = ""

    for nt in seq:
        current_index = nucleotide_to_index[nt]
        probabilities = calculate_p_matrix(t=t, mu=mu)[current_index]
        evolved_nt_index = np.random.choice([0, 1, 2, 3], p=probabilities)
        evolved_sequence += index_to_nucleotide[evolved_nt_index]

    return evolved_sequence


def calculate_steps_for_mutations(mutation_percentage, total_nucleotides, mu):
    target_mutations = (mutation_percentage / 100) * total_nucleotides
    expected_mutations_per_step = total_nucleotides * mu
    return target_mutations / expected_mutations_per_step


# Calculate the number of the percentage point change from the tree and then use it to find our the number of expected
# steps to use in the simulation.
TREE_LEAVES_JC_STEPS = {}
for leaf in original_tree.get_terminals():
    distance_to_root = original_tree.distance(leaf)
    print(f"Distance from {leaf.name} to the root: {distance_to_root}")
    steps_needed = int(
        np.round(calculate_steps_for_mutations(distance_to_root, total_nucleotides=len(back_translated_sequence),
                                               mu=1.0 / 5000)))
    print(f"Steps needed to mutate original sequence to {leaf.name}: {steps_needed}")
    TREE_LEAVES_JC_STEPS.update({leaf.name: steps_needed})

RANDOMLY_GENERATED_TREE_SEQUENCES = {}
# Now just run the Markov JC with the appropriate steps.
for leaf, timesteps in TREE_LEAVES_JC_STEPS.items():
    simulated_mutation_sequence = markov_jukes_cantor(seq=back_translated_sequence, t=timesteps)
    RANDOMLY_GENERATED_TREE_SEQUENCES.update({leaf: Seq(simulated_mutation_sequence)})


# Sanity Check --- Calculate the average percent point change in the simulated sequences
# I know that it is not supposed to be perfectly accurate due to the stochastic nature of the model.
def percent_hamming_distance(original, modified):
    distance = sum(ch1 != ch2 for ch1, ch2 in zip(original, modified))
    return 100 * (distance / len(original))


for leaf, seq in RANDOMLY_GENERATED_TREE_SEQUENCES.items():
    print(
        f"Percentage of point mutations between Original and {leaf} : {percent_hamming_distance(back_translated_sequence, seq)}")

# Its very close so I assume correctness


# Step 4:
# Translate the generated sequences into proteins and save the generated sequences in a fasta file (proteins.fa).

protein_sequences = [seq.translate(table='Standard') for _, seq in RANDOMLY_GENERATED_TREE_SEQUENCES.items()]
protein_records = [SeqRecord(protein_seq, id=f"{leaf}", description="") for leaf, protein_seq in
                   zip(RANDOMLY_GENERATED_TREE_SEQUENCES.keys(), protein_sequences)]

with open("../data/proteins.fa", "w") as output_handle:
    SeqIO.write(protein_records, output_handle, "fasta")

# Step 5:
# Create a multiple alignment based on these protein sequences and the BLOSUM matrix using the clustalw program.

# Note: I have a Windows machine, so I downloaded clustalw2.exe and connected it to the BioPython Applications
# I do not know if this code works out-of-the-box if it is run to another machine.
# I guess since there is a ClustalW GUI available online anyone can verify the alignment given the proteins file.
clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="../data/proteins.fa", matrix="BLOSUM", output='NEXUS')
subprocess.call(str(clustalw_cline), shell=True)


# Step 6:
# Based on the multiple alignment, generate phylogenetic trees using the upgma and nj (neighbor-joining) methods.
alignment = AlignIO.read('../data/proteins.nxs', 'nexus')

# Protein calculator distance with ‘blosum62’ model
def get_distance_matrix(aln):
    calculator = DistanceCalculator('blosum62')
    return calculator.get_distance(aln)


dist_mat = get_distance_matrix(alignment)
constructor = DistanceTreeConstructor()

# UPGMA METHOD
upgmatree = constructor.upgma(dist_mat)

# NJ METHOD
njtree = constructor.nj(dist_mat)

# Let's draw them
Phylo.draw(upgmatree, branch_labels=lambda c: round(c.branch_length, 2))
Phylo.draw(njtree, branch_labels=lambda c: round(c.branch_length, 2))

# And finally save them
Phylo.write(upgmatree, "../data/upgmatree.nxs", "nexus")
Phylo.write(njtree, "../data/njtree.nxs", "nexus")