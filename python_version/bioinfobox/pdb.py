from Bio.PDB import PDBParser, NeighborSearch, Structure, Atom
from itertools import combinations
from numpy import float32


def load_pdb(path: str, name: str) -> Structure.Structure:
    """
    Parses a PDB file.

    :param path: The path to the file.
    :param name: User provided name for the structure.
    :return: The structure encoded by the input file.
    """
    parser = PDBParser()
    structure = parser.get_structure(name, path)
    return structure


def pdb_get_models(pdb_structure: Structure.Structure):
    """
    Returns all models.

    :param pdb_structure: A structure parsed from the PDB format.
    :return: A class representing the models of the structure.
    """
    return pdb_structure.get_models()


def pdb_get_chains(pdb_structure: Structure.Structure):
    """
    Returns all chains.

    :param pdb_structure: A structure parsed from the PDB format.
    :return: A class representing the chains of the structure.
    """
    return pdb_structure.get_chains()


def pdb_get_residues(pdb_structure: Structure.Structure):
    """
    Returns all of the residues.

    :param pdb_structure: A structure parsed from the PDB format.
    :return: A class representing the residues of the structure.
    """
    return pdb_structure.get_residues()


def pdb_get_atoms(pdb_structure: Structure.Structure):
    """
    Returns all the atoms of the structure.

    :param pdb_structure: A structure parsed from the PDB format.
    :return: A class representing the atoms of the structure.
    """
    return pdb_structure.get_atoms()


def pdb_info(pdb_structure: Structure.Structure) -> dict:
    """
    Creates a dictionary of lengths of models, chains, etc.

    :param pdb_structure: A structure parsed from the PDB format.
    :return: A dictionary containing lengths of models, chains, etc. of the structure.
    """
    info_dict = {'models': len(list(pdb_get_models(pdb_structure))),
                 'chains': len(list(pdb_get_chains(pdb_structure))),
                 'residues': len(list(pdb_get_residues(pdb_structure))),
                 'atoms': len(list(pdb_get_atoms(pdb_structure)))}
    return info_dict


def pdb_width(pdb_structure: Structure.Structure) -> float32:
    """
    Returns the maximum distance between two atoms in the structure.

    :param pdb_structure: A structure parsed from the PDB format.
    :return: The maximum distance between atoms in the structure.
    """
    atoms = pdb_get_atoms(pdb_structure)
    atomic_pairs = combinations(atoms, 2)
    width = max([b - a for (a, b) in atomic_pairs])
    return width


def pdb_get_neighbor_atoms(hetatm: Atom.Atom, distance: float, pdb_structure: Structure.Structure) -> list:
    """
    Returns a list of atoms in a defined distance from a given atom.

    :param hetatm: An atom from which the distance is measured.
    :param distance: A radius we are considering.
    :param pdb_structure: The structure in which we are searching.
    :return: A list of atoms.
    """
    atoms = list(pdb_get_atoms(pdb_structure))
    neighborhood_search = NeighborSearch(atoms)
    return neighborhood_search.search(hetatm.get_coord(), distance, level="A")


def pdb_get_neighbor_residues(hetatm: Atom.Atom, distance: float, pdb_structure: Structure.Structure) -> list:
    """
    Returns a list of residues in a defined distance from a given atom.

    :param hetatm: An atom from which the distance is measured.
    :param distance: A radius we are considering.
    :param pdb_structure: The structure in which we are searching.
    :return: A list of residues.
    """
    atoms = list(pdb_get_atoms(pdb_structure))
    neighborhood_search = NeighborSearch(atoms)
    return neighborhood_search.search(hetatm.get_coord(), distance, level="R")
