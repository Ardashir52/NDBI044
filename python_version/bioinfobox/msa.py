from Bio import SeqIO, Seq
from itertools import combinations
from Bio.Align import substitution_matrices


def load_clustal(path: str) -> list:
    """
    Parses a clustal file.

    :param path: The path to the file to parse.
    :return: A list of aligned sequences.
    """
    return list(SeqIO.parse(path, "clustal"))


def clustal_retrieve_position(clustal_list: list, position: int) -> Seq.Seq:
    """
    Returns a sequence at a given position.

    :param clustal_list: The parsed MSA as a list.
    :param position: Position of the desired sequence.
    :return: The chosen sequence.
    """
    return clustal_list[position].seq


def clustal_retrieve_id(clustal_list: list, clustal_id: str) -> Seq.Seq:
    """
    Returns a sequence determined by the given id.

    :param clustal_list: The parsed MSA as a list.
    :param clustal_id: The id of desired sequence.
    :return: The chosen sequence.
    """
    for clustal_record in clustal_list:
        if clustal_record.id == clustal_id:
            return clustal_record.seq


def clustal_column(clustal_list: list, position: int) -> list:
    """
    Returns a chosen column of the MSA.

    :param clustal_list: The parsed MSA as a list.
    :param position: The number of the chosen column.
    :return: The column at a given position as a list.
    """
    column = []
    for cl_record in clustal_list:
        column.append(cl_record.seq[position])
    return column


def sum_of_column(column: list, sub_mx: substitution_matrices.Array) -> int:
    """
    Computes the sum of pairs score of a column.

    :param column: A column from a MSA.
    :param sub_mx: Chosen substitution matrix.
    :return: The sum of pairs score of the column.
    """
    combos = list(combinations(column, 2))
    column_sum = 0
    alphabet = sub_mx.alphabet
    for pair in combos:
        if pair[0] not in alphabet and pair[0] == pair[1]:
            column_sum += 1
        elif pair[0] not in alphabet or pair[1] not in alphabet:
            column_sum += -1
        else:
            column_sum += sub_mx[pair[0], pair[1]]
    return column_sum


def msa_score(clustal_list: list, sub_mx: substitution_matrices.Array) -> int:
    """
    Computes the sum of pairs score over all the columns of the MSA.

    :param clustal_list: The parsed MSA as a list.
    :param sub_mx: Chosen substitution matrix.
    :return: The sum of pairs score of the MSA.
    """
    msa_sum = 0
    for position in range(len(clustal_list)):
        column = clustal_column(clustal_list, position)
        msa_sum += sum_of_column(column, sub_mx)
    return msa_sum
