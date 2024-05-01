import os
import math
import numpy as np
from hmmlearn import hmm
from Bio import SeqIO
from itertools import product
import pickle

### Get the training data (cpg.fa file)
TRAIN_DATA = [str(f.seq) for f in SeqIO.parse('../data/cpg.fa', 'fasta')][0]

# We know that we must use a HMM with 2 States  Cpg(+) (CpG island) / CpG(-) (no CpG island)
# Regarding the choice of the emissions we see that choosing a single nucleotide does not
# encapsulate the problem well. So we will use dinucleotide sequences as our emissions instead.
nucleotides = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
dinucleotides = {k: v for v, k in enumerate([''.join(pair) for pair in product(nucleotides.keys(), repeat=2)])}


def encode_sequence(seq):
    """
    Encodes a Sequence into a dinucleotide
    :param seq: A BioPython Seq object
    :return: A train-ready numpy array in (len(seq), 1) shape as the CategoricalHMM expects
    """
    encoded = np.array([dinucleotides[seq[i - 1:i + 1]] for i in range(1, len(seq), 2)]).reshape(-1, 1)
    return encoded


def save_model(model, path='../data/cat_hmm.pt'):
    """
    Saves a Categorical HMM model by saving its params
    """
    params = {}
    params.update({'startprob_': model.startprob_})
    params.update({'transmat_': model.transmat_})
    params.update({'emissionprob_': model.emissionprob_})
    with open(path, 'wb') as fout:
        pickle.dump(params, fout)
    return


def load_model(model, path='../data/cat_hmm.pt'):
    """
    Loads a Categorical HMM model from a pickle object
    """
    with open(path, 'rb') as fin:
        params = pickle.load(fin)
    model.startprob_ = params['startprob_']
    model.transmat_ = params['transmat_']
    model.emissionprob_ = params['emissionprob_']
    return


train_data_enc = encode_sequence(TRAIN_DATA)
hmm_model = hmm.CategoricalHMM(n_components=2, n_iter=100)
if os.path.exists('../data/cat_hmm.pt'):
    load_model(hmm_model, '../data/cat_hmm.pt')
else:
    hmm_model.fit(train_data_enc)

#### Print the results
print(f"Final Transition Matrix: {hmm_model.transmat_}")
print(f"Final Emission Matrix: {hmm_model.emissionprob_}")





### Get the test data (cpg_test.fa file)
TEST_DATA = [str(f.seq) for f in SeqIO.parse('../data/cpg_test.fa', 'fasta')]


### There are equivalent ways of calculating posterior probabilities ###
### We can set the state to 0 and just calculate a chain of products from the emission table
### This is due to the fact that we assume the model stays permanently in state 0 and doesnt transition
### Then we can to the same for class 1.
### I tried to do that but i get underflow errors, so the best way is to calculate log-probs
### This is equivalent to the methods below which return pretty much the same results


def posterior_prob_by_most_probable_states_hard(seq, model):
    """
    Returns the probability of a sequence coming from a state 0 (CpG Island) vs state 1 (no Island)
    For this calculation this function uses the rates of most probable classes as returned by the
    forward-backward algorithm.
    """
    hard_class_assignment = model.predict(seq)
    prob_class_0 = round(hard_class_assignment.sum() / len(hard_class_assignment), 2)
    prob_class_1 = 1 - prob_class_0
    return prob_class_0, prob_class_1


def posterior_prob_by_most_probable_states_soft(seq, model):
    """
    Returns the probability of a sequence coming from a state 0 (CpG Island) vs state 1 (no Island)
    For this calculation this function averages the scores over the per-state probabilities the model gets
    from the forward-backward algorithm.
    """
    soft_class_assignment = model.predict_proba(seq)
    prob_class_1 = round(soft_class_assignment.mean(axis=0)[1], 2)
    prob_class_0 = 1 - prob_class_1
    return prob_class_0, prob_class_1


### Sanity Check ###
encoded_seq = encode_sequence('ATATAATATAATATAATATAATATAATATA')
soft_class_0_prob, soft_class_1_prob = posterior_prob_by_most_probable_states_soft(encoded_seq, hmm_model)
print(soft_class_0_prob, soft_class_1_prob) # Returns 0.0 / 1.0

encoded_seq = encode_sequence('CGCGCGCGCGCGCGCGCGCGCGCGCGCGCG')
soft_class_0_prob, soft_class_1_prob = posterior_prob_by_most_probable_states_soft(encoded_seq, hmm_model)
print(soft_class_0_prob, soft_class_1_prob) # Returns 0.9 / 0.1


with open("../data/results.txt", "w") as output_file:
    for seq in TEST_DATA:
        encoded_seq = encode_sequence(seq)
        soft_class_0_prob, soft_class_1_prob = posterior_prob_by_most_probable_states_soft(encoded_seq, hmm_model)
        # Write the result to a file
        output_file.write(f" Non-CpG Island / CpG Island: {soft_class_0_prob, soft_class_1_prob}\n")

print("Posterior probabilities have been written to results.txt.")
