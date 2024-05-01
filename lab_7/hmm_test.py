from hmmlearn import hmm
import numpy as np
from Bio import SeqIO


with open('../data/cpg.fa') as f:
    for r in SeqIO.parse(f, 'fasta'):
        seq = str(r.seq)
print(seq[:20])

nucleotides = ['A', 'T', 'C', 'G']
pairs_nt = []
for nt1 in nucleotides:
    for nt2 in nucleotides:
        pairs_nt.append(nt1 + nt2)

Y = [pairs_nt.index(seq[i - 1:i + 1]) for i in range(1, len(seq), 2)]
Y = np.array(Y)

print(Y[:20])

model2 = hmm.CategoricalHMM(n_components=2, n_iter=100)
model2.fit(Y.reshape(-1, 1))

print("Macierz przejść:")
print(model2.transmat_)

print("Macierz emisji:")
print(model2.emissionprob_)

T_y = model2.transmat_
print(T_y)
E_y = 1 / T_y[0][1]
print(E_y)





