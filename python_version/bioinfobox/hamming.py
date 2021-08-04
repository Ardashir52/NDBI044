from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def hamming(seq1, seq2) -> int:
    """
    Computes the Hamming distance of two sequences, that is the number
    of position where two sequences of equal length differ.

    :param seq1: First sequence.
    :param seq2: Second sequence.
    :return: The Hamming distance of the two sequences.
    """
    if type(seq1) is SeqRecord:
        return hamming(seq1.seq, seq2)
    elif type(seq2) is SeqRecord:
        return hamming(seq1, seq2.seq)
    elif (type(seq1) is str or type(seq1) is Seq) and (type(seq2) is Seq or type(seq2) is str):
        if len(seq1) != len(seq2):
            raise ValueError('The sequences are of different lengths!')
        else:
            distance = 0
            for i in range(len(seq1)):
                if seq1[i] != seq2[i]:
                    distance += 1
        return distance
    else:
        raise TypeError('Wrong type.')
