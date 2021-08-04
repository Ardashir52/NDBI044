from Bio.Align import substitution_matrices


def load_matrix(name: str) -> substitution_matrices.Array:
    """
    Loads a substitution matrix from the ones built in the Biopython library.

    :param name: The name of the matrix to load.
    :return: The chosen substitution matrix.
    """
    return substitution_matrices.load(name)


def matrix_offer():
    """
    Provides a list of available substitution matrices.

    :return: Prints a numbered list of matrices.
    """
    mxs = substitution_matrices.load()
    for elem in range(len(mxs)):
        print(str((elem + 1)) + " " + mxs[elem])


def matrix_choice(num: int) -> str:
    """
    Returns a name of the matrix based on its number in the list
    from the matrix_offer function.

    :param num: The position of the matrix in the list from the matrix_offer function.
    :return: Name of the matrix at a given position.
    """
    mxs = substitution_matrices.load()
    return mxs[num - 1]
