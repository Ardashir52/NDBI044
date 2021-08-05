from statistics import mean
from Bio.PDB import Atom, Structure
from bioinfobox.pdb import pdb_get_neighbor_atoms
from bioinfobox.conservation import conservation_at_position


def structure_sequence_equivalence(hetatm: Atom.Atom, clustal_msa: list,
                                   pdb_structure: Structure.Structure, offset=0) -> (float, float):
    """
    Checks whether active sites are indeed more conserved.

    :param hetatm: A ligand determining the active site in the atom format.
    :param clustal_msa: A parsed MSA in a form of a list.
    :param pdb_structure: A parsed PDB structure in which we search.
    :param offset: Defines a shift in the sequence; 0 by default.
    :return: A tuple of 1. conservation score of the active site, 2. conservation score of the whole alignment.
    """
    neighbors = pdb_get_neighbor_atoms(hetatm, 3, pdb_structure)
    positions = [atom.full_id[3][1] for atom in neighbors]
    cons_scores_at_positions = [conservation_at_position(clustal_msa, pos + offset) for pos in positions]
    seq_len = len(clustal_msa[0].seq)
    all_scores = [conservation_at_position(clustal_msa, i) for i in range(seq_len)]
    return mean(cons_scores_at_positions), mean(all_scores)
