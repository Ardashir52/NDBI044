from Bio import SeqIO, Seq, SeqRecord


def fasta_get_seq_and_desc(record: SeqRecord.SeqRecord) -> (Seq.Seq, str):
    """
    Returns the sequence and its description.

    :param record: Data of a single molecule in the SeqRecord type.
    :return: Tuple of the sequence and its description.
    """
    return record.seq, record.description


def fasta_seq_length(record: SeqRecord.SeqRecord) -> int:
    """
    Returns the length of the sequence in the SeqRecord provided.

    :param record: Data of a single molecule in the SeqRecord type.
    :return: Length of the sequence.
    """
    return len(record.seq)


def fasta_seq_subsequence(record: SeqRecord.SeqRecord, start: int, end: int) -> Seq.Seq:
    """
    Return a subsequence of the input sequence based on start
    and end parameters.

    :param record: Data of a single molecule in the SeqRecord type.
    :param start: First position of the subsequence.
    :param end: First position after the subsequence.
    :return: Subsequence of the input sequence
    """
    return record.seq[start:end]


def load_fasta(path: str) -> dict:
    """
    Loads provided fasta file to a dictionary of sequence records.
    Presumed is the uniqueness of their keys.

    :param path: Path to the FASTA file.
    :return: Dictionary of the sequence records.
    """
    return SeqIO.to_dict(SeqIO.parse(path, "fasta"))
