import itertools
import sys
from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.Align import substitution_matrices

BLOSUM80 = substitution_matrices.load("BLOSUM80")

### Read the DNA Sequences
histones_records = [f.seq for f in list(SeqIO.parse("../data/histones.fa", "fasta"))]
bzips_records = [f.seq for f in list(SeqIO.parse("../data/bzips.fa", "fasta"))]


def align_pair(s1, s2, mode, protein):
    aligner = Align.PairwiseAligner()

    try:
        aligner.mode = mode
    except:
        raise Exception

    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    if protein:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")

    return aligner.align(s1, s2)


def get_combinations_of_2(iterable):
    return list(itertools.combinations(iterable, 2))


def calc_alignment_seq2seq(s1, s2, matrix=None, translate=False):
    # Global alignment 1 for match -1 for mismatch, -1 open break penalty, -0.5 expand break penalty#
    if translate:
        s1 = standard_translation(s1)
        s2 = standard_translation(s2)
    if matrix is not None:
        g_score = pairwise2.align.globalds(s1, s2, matrix, -1, -0.5)[0].score
        l_score = pairwise2.align.globalds(s1, s2, matrix, -1, -0.5)[0].score
    else:
        g_score = pairwise2.align.globalms(s1, s2, 1, -1, -1, -0.5)[0].score
        g_score2 = align_pair(s1, s2, mode='global', protein=False).score
        assert g_score == g_score2
        l_score = pairwise2.align.localms(s1, s2, 1, -1, -1, -0.5)[0].score
        l_score2 = align_pair(s1, s2, mode='local', protein=False).score
        assert l_score == l_score2
    return g_score, l_score


def calc_alignment_in_group(group, matrix=None, translate=False):
    g_scores = 0.
    l_scores = 0.
    combs = get_combinations_of_2(group)
    for comb in combs:
        g, l = calc_alignment_seq2seq(comb[0], comb[1], matrix=matrix, translate=translate)
        g_scores += g
        l_scores += l

    g_scores /= len(combs)
    l_scores /= len(combs)
    return g_scores, l_scores



def calc_alignment_between_groups(group_1, group_2, matrix=None, translate=False):
    g_scores = 0.
    l_scores = 0.
    for g1 in group_1:
        for g2 in group_2:
            g, l = calc_alignment_seq2seq(g1, g2, matrix=matrix, translate=translate)
            g_scores += g
            l_scores += l

    g_scores /= (len(group_1) + len(group_2))
    l_scores /= (len(group_1) + len(group_2))
    return g_scores, l_scores


def standard_translation(sequence):
    return sequence.translate(table='Standard')


print("Non-Translated")
#histones_g, histones_l = calc_alignment_in_group(histones_records)
# histones_g = 112.80
# histones_l = 118.89
#print(f"Histones Score: Global {hist_g} / Local {hist_l}")

# # bzip_g, bzip_l = calc_alignment_in_group(bzips_records)
# bzip_g = 101.21
# bzip_l = 176.60
# print(f"Bzip Score: Global {bzip_g} / Local {bzip_l}")
#
# # between_g, between_l = calc_alignment_between_groups(histones_records, bzips_records)
# between_g = -515.78
# between_l = 424.66
# print(f"Between Group Score: Global {between_g} / Local {between_l}")
#
# print("Translated")
# histones_g, histones_l = calc_alignment_in_group(histones_records, translate=True)
# # histones_g = 112.80
# # histones_l = 118.89
# print(f"Histones Score: Global {histones_g} / Local {histones_l}")
#
# bzip_g, bzip_l = calc_alignment_in_group(bzips_records, True)
# # bzip_g = 101.21
# # bzip_l = 176.60
# print(f"Bzip Score: Global {bzip_g} / Local {bzip_l}")
#
# between_g, between_l = calc_alignment_between_groups(histones_records, bzips_records, True)
# # between_g = -515.78
# # between_l = 424.66
# print(f"Between Group Score: Global {between_g} / Local {between_l}")
