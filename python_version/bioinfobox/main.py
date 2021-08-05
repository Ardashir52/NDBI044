import argparse
import sys

from bioinfobox.ed_align import edit_distance
from bioinfobox.hamming import hamming
from bioinfobox.fasta import *
from bioinfobox.msa import load_clustal, clustal_retrieve_position, clustal_retrieve_id, clustal_column, msa_score
from bioinfobox.conservation import conservation_at_position, n_most_conserved
from bioinfobox.pdb import load_pdb, pdb_width
from bioinfobox.matrix import load_matrix, matrix_offer


def frontend():
    """
    Provides a command line interface.
    """
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    parser_edit = subparser.add_parser('edit', help="Computes edit distance and finds all alignments.")
    parser_edit.set_defaults(action='edit')
    parser_edit.add_argument('seq1')
    parser_edit.add_argument('seq2')
    parser_edit.add_argument('--align', action='store_true')

    parser_hamming = subparser.add_parser('hamming', help="Computes hamming distance of two sequences.")
    parser_hamming.set_defaults(action='hamming')
    parser_hamming.add_argument('seq1')
    parser_hamming.add_argument('seq2')

    parser_fasta = subparser.add_parser('fasta', help="Parses a fasta file and provides basic info.")
    parser_fasta.set_defaults(action='fasta')
    parser_fasta.add_argument('path')
    parser_fasta.add_argument('--sequence', dest='seq', type=str)
    parser_fasta.add_argument('--description', dest='desc', type=str)
    parser_fasta.add_argument('--len', dest='len', type=str)
    parser_fasta.add_argument('--subsequence', dest='subseq', nargs=3)

    parser_msa = subparser.add_parser('msa', help="Parses a clustal file and provides simple operations thereover.")
    parser_msa.set_defaults(action='msa')
    parser_msa.add_argument('path')
    parser_msa.add_argument('--position', dest='pos', type=int)
    parser_msa.add_argument('--id', dest='id', type=str)
    parser_msa.add_argument('--column', dest='column', type=int)
    parser_msa.add_argument('--score', dest='score', type=str)

    parser_conservation = subparser.add_parser('conservation', help="Determines the conservation of a position.")
    parser_conservation.set_defaults(action='conservation')
    parser_conservation.add_argument('path')
    parser_conservation.add_argument('--position', dest='position', type=int)
    parser_conservation.add_argument('--best', dest='n', type=int)

    parser_pdb = subparser.add_parser('pdb', help="Parser pdb file and provides additional info.")
    parser_pdb.set_defaults(action='pdb')
    parser_pdb.add_argument('path')
    # parser_pdb.add_argument('--models', action='store_true')
    # parser_pdb.add_argument('--chains', action='store_true')
    # parser_pdb.add_argument('--residues', action='store_true')
    # parser_pdb.add_argument('--atoms', action='store_true')
    # parser_pdb.add_argument('--width', action='store_true')
    # parser_pdb.add_argument('--neigh_atom', dest='neigh_atom', nargs=6)
    # parser_pdb.add_argument('--neigh_res', dest='neigh_res', nargs=6)

    # parser_str_seq = subparser.add_parser('str_seq', help="Provides comparison of conservation of an active site to "
    #                                                      "the whole sequence.")
    # parser_str_seq.set_defaults(action='str_seq')
    # parser_str_seq.add_argument('pdb_path')
    # parser_str_seq.add_argument('clustal_path')
    # parser_str_seq.add_argument('-m', dest='model', required=True)
    # parser_str_seq.add_argument('-c', dest='chain', required=True)
    # parser_str_seq.add_argument('-r', dest='residue', required=True, nargs=3)
    # parser_str_seq.add_argument('-a', dest='atom', required=True)
    # parser_str_seq.add_argument('--offset', dest='offset', type=int)

    parser_matrix = subparser.add_parser('matrix', help="Provides help for picking a substitution matrix.")
    parser_matrix.set_defaults(action='matrices')

    if len(sys.argv) < 2:
        parser.print_help()
        parser.print_usage()
    else:
        choices = parser.parse_args()
        if choices.action == 'edit':
            num, align_list = edit_distance(choices.seq1, choices.seq2)
            print("The edit distance is: " + str(num))
            if choices.align:
                for el in align_list:
                    print()
                    print(el[0])
                    print(el[1])
        elif choices.action == 'hamming':
            print("The Hamming distance of the sequences is: " + str(hamming(choices.seq1, choices.seq2)))
        elif choices.action == 'fasta':
            fasta_dict = load_fasta(choices.path)
            if choices.desc is not None:
                _, desc = fasta_get_seq_and_desc(fasta_dict[choices.desc])
                print(desc)
            if choices.seq is not None:
                seq, _ = fasta_get_seq_and_desc(fasta_dict[choices.seq])
                print(seq)
            if choices.len is not None:
                length = fasta_seq_length(fasta_dict[choices.len])
                print("The length of the sequence is: " + str(length))
            if choices.subseq is not None:
                print(choices.subseq)
                subseq = fasta_seq_subsequence(fasta_dict[choices.subseq[0]], int(choices.subseq[1]),
                                               int(choices.subseq[2]))
                print(subseq)
        elif choices.action == 'msa':
            clustal_list = load_clustal(choices.path)
            if choices.pos is not None:
                print(clustal_retrieve_position(clustal_list, choices.pos))
            if choices.id is not None:
                print(clustal_retrieve_id(clustal_list, choices.id))
            if choices.column is not None:
                col = clustal_column(clustal_list, choices.column)
                print(col)
            if choices.score is not None:
                mx = load_matrix(choices.score)
                print(msa_score(clustal_list, mx))
        elif choices.action == 'conservation':
            clustal_list = load_clustal(choices.path)
            if choices.position is not None:
                cons = conservation_at_position(clustal_list, choices.position)
                print("The conservation score of the column is: " + str(cons))
            if choices.n is not None:
                n_cons = n_most_conserved(clustal_list, choices.n)
                print(n_cons)
        elif choices.action == 'pdb':
            pdb_structure = load_pdb(path=choices.path, name='str1')
            width = pdb_width(pdb_structure)
            print("The width of the structure is: " + str(width))
        elif choices.action == 'matrices':
            matrix_offer()


if __name__ == '__main__':
    frontend()
