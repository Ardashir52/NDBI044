from statistics import mode
from bioinfobox.msa import clustal_column


def is_residue(symbol: str) -> bool:
    """
    Checks if an input is a valid symbol for a residue.

    :param symbol: A symbol to check.
    :return: A judgement on whether the given symbol is a residue.
    """
    return symbol in "ARNDCQEGHILKMFPSTWYV"


def conservation_value(column: list) -> float:
    """
    Computes the conservation value of a given column.

    :param column: A column to evaluate.
    :return: The conservation value of the column.
    """
    col = list(filter(is_residue, column))
    return col.count(mode(col)) / len(column)


def conservation_at_position(clustal_list: list, position: int) -> float:
    """
    Returns a conservation score of a column at given position.

    :param clustal_list: A parsed MSA in the form of a list.
    :param position: The position of the column to evaluate.
    :return: A conservation score of the chosen column.
    """
    column = clustal_column(clustal_list, position)
    return conservation_value(column)


def n_most_conserved(clustal_list: list, n: int) -> list:
    """
    Returns a list of the n most conserved positions for a given n and a MSA.

    :param clustal_list: A parsed MSA in the form of a list.
    :param n: The number of the most conserved positions we want to obtain.
    :return: A list of n most conserved position.
    """
    seq_len = len(clustal_list[0].seq)
    cons_list = [(i + 1, conservation_at_position(clustal_list, i)) for i in range(seq_len)]
    return sorted(cons_list, key=lambda pair: pair[1], reverse=True)[:n]
