from nose.tools import eq_
from pathlib import Path
from os.path import join

from bioinfobox.hamming import hamming
from bioinfobox.fasta import load_fasta, fasta_seq_subsequence
from bioinfobox.ed_align import edit_distance
from bioinfobox.msa import load_clustal, clustal_retrieve_id, clustal_retrieve_position, msa_score
from bioinfobox.conservation import n_most_conserved
from bioinfobox.pdb import load_pdb, pdb_info, pdb_get_neighbor_atoms
from bioinfobox.str_seq import structure_sequence_equivalence
from bioinfobox.matrix import load_matrix

PATH = Path(__file__).parent.absolute()


def test_fasta_subseq():
    fasta_list = load_fasta(join(PATH, "fastaseqtest"))
    rec = fasta_list['MCHU']
    eq_(fasta_seq_subsequence(rec, 0, 5), 'MADQL')


def test_hamming():
    fasta_list = load_fasta(join(PATH, "fastaseqtest"))
    rec = fasta_list['MCHU']
    substr = fasta_seq_subsequence(rec, 0, 5)
    eq_(hamming(substr, 'MADOL'), 1)


def test_ed():
    num, align_list = edit_distance('Sunday', 'Saturday')
    eq_(num, 3)
    eq_(len(align_list), 3)


def test_load_clustal():
    rec = load_clustal(join(PATH, "p53_mafft_clustal.txt"))
    eq_(clustal_retrieve_id(rec, "UniRef90_G3WS63"), clustal_retrieve_position(rec, 54))


def test_msa():
    rec = load_clustal(join(PATH, "p53_mafft_clustal.txt"))
    mx = load_matrix('BLOSUM62')
    eq_(msa_score(rec, mx), 83756.0)


def test_conservation_score():
    rec = load_clustal(join(PATH, "p53_mafft_clustal.txt"))
    eq_(n_most_conserved(rec, 5)[3][1], 1.0)


def test_pdb_info():
    pdb_str = load_pdb(join(PATH, 'pdb1lb5.ent.txt'), '1lb5')
    eq_(pdb_info(pdb_str)['residues'], 242)


def test_pdb_neighbors():
    pdb_str = load_pdb(join(PATH, 'pdb1lb5.ent.txt'), '1lb5')
    neighborhood = pdb_get_neighbor_atoms(pdb_str[0]['B'][('W', 821, ' ')]['O'], 5, pdb_str)
    eq_(len(neighborhood), 6)


def test_sequence_structure():
    pdb_str = load_pdb(join(PATH, 'pdb3kz8.ent.txt'), '3kz8')
    seq = load_clustal(join(PATH, "p53_mafft_clustal (kopie).txt"))
    a, b = structure_sequence_equivalence(pdb_str[0]['B'][('H_ZN', 1, ' ')]['ZN'], seq, pdb_str, 80)
    eq_(a > b, True)
