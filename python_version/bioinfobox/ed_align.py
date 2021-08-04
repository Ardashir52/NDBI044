from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def ed_backtrack(position: (int, int), partial_word1, partial_word2, seq_list: list, table: [[int]], seq1, seq2):
    """
    A function that builds the alignments for the edit_distance function by backtracking over
    the filled table.

    :param position: Current position in the table, used as an end condition for recursion.
    :param partial_word1: Progressively built alignment.
    :param partial_word2: Progressively built alignment.
    :param seq_list: A list of all alignments so far found.
    :param table: The alignment table precomputed in the edit_distance function.
    :param seq1: The first sequence compared.
    :param seq2: The second sequence compared.
    :return: Nothing.
    """
    if position == (0, 0):
        w1 = partial_word1[::-1]
        w2 = partial_word2[::-1]
        seq_list.append((w1, w2))
    else:
        coord1 = position[0]
        coord2 = position[1]
        minimum = min(table[coord1 - 1][coord2 - 1], table[coord1 - 1][coord2], table[coord1][coord2 - 1])
        if table[coord1 - 1][coord2 - 1] == minimum:
            new_word1 = partial_word1 + seq1[coord1 - 1]
            new_word2 = partial_word2 + seq2[coord2 - 1]
            ed_backtrack((coord1 - 1, coord2 - 1), new_word1, new_word2, seq_list, table, seq1, seq2)
        if table[coord1 - 1][coord2] == minimum:
            new_word1 = partial_word1 + seq1[coord1 - 1]
            new_word2 = partial_word2 + "–"
            ed_backtrack((coord1 - 1, coord2), new_word1, new_word2, seq_list, table, seq1, seq2)
        if table[coord1][coord2 - 1] == minimum:
            new_word1 = partial_word1 + "–"
            new_word2 = partial_word2 + seq2[coord2 - 1]
            ed_backtrack((coord1, coord2 - 1), new_word1, new_word2, seq_list, table, seq1, seq2)


def edit_distance(seq1, seq2) -> (int, list):
    """
    Function for measuring edit distance and finding all corresponding alignments.

    :param seq1: A sequence to compare.
    :param seq2: A sequence to compare.
    :return: Value of the distance and a list of all alignments (as tuples of sequences).
    """
    if type(seq1) is SeqRecord:
        return edit_distance(seq1.seq, seq2)
    elif type(seq2) is SeqRecord:
        return edit_distance(seq1, seq2.seq)
    elif (type(seq1) is str or type(seq1) is Seq) and (type(seq2) is Seq or type(seq2) is str):
        table = [[0 for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
        for i in range(len(table)):
            table[i][0] = i
        for j in range(len(table[0])):
            table[0][j] = j
        for j in range(1, len(table[0])):
            for i in range(1, len(table)):
                if seq1[i - 1] == seq2[j - 1]:
                    cost = 0
                else:
                    cost = 1
                table[i][j] = min(table[i - 1][j] + 1, table[i][j - 1] + 1, table[i - 1][j - 1] + cost)
        alignment_list = []
        ed_backtrack((len(seq1), len(seq2)), "", "", alignment_list, table, seq1, seq2)
        return table[len(seq1)][len(seq2)], alignment_list
    else:
        raise TypeError('Wrong type.')
