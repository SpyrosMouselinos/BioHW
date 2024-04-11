from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator

alignment_PAH = AlignIO.read('../data/Human_PAH_paralogues.nex', 'nexus')
alignment_H2BFS = AlignIO.read('../data/Human_H2BFS_paralogues.nex', 'nexus')
alignment_PAH_o = AlignIO.read('../data/Human_PAH_orthologues.nex', 'nexus')

print(alignment_PAH)
print(alignment_H2BFS)


# Protein calculator distance with ‘blosum62’ model
def get_distance_matrix(aln):
    calculator = DistanceCalculator('blosum62')
    return calculator.get_distance(aln)


dist_m_PAH = get_distance_matrix(alignment_PAH)
dist_m_H2BFS = get_distance_matrix(alignment_H2BFS)
# dist_m_PAH_o = get_distance_matrix(alignment_PAH_o)


constructor = DistanceTreeConstructor()

upgmatree_PAH = constructor.upgma(dist_m_PAH)
njtree_PAH = constructor.nj(dist_m_PAH)

upgmatree_H2BFS = constructor.upgma(dist_m_H2BFS)
njtree_H2BFS = constructor.nj(dist_m_H2BFS)

Phylo.draw(upgmatree_PAH)
Phylo.draw(njtree_PAH)

Phylo.draw(upgmatree_H2BFS, branch_labels=lambda c: c.branch_length)
Phylo.draw(njtree_H2BFS, branch_labels=lambda c: round(c.branch_length, 3))
