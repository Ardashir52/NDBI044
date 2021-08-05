# NDBI044

##Bioinformatics toolbox

This package implements some basic bioinformatic tools for simple tasks.

###Installation

####For installation using python virtual environment follow the following script:
```console
$ git clone https://github.com/Ardashir52/NDBI044.git
$ cd python_version
$ python3 -m venv venv
$ source venv/bin/activate
$ python3 setup.py install
```

####For a global install follow just these steps:
```console
$ git clone https://github.com/Ardashir52/NDBI044.git
$ cd python_version
$ python3 setup.py install
``` 

##Library functionality

###Edit distance

This function allows to compute the edit distance of two sequences and to retreive
all of the corresponding alignments.

Example:

```
>>> from bioinfobox.ed_align import edit_distance
>>> distance, list_of_aligns = edit_distance('Sunday', 'Saturday')
>>> print(distance, list_of_aligns)
3, [('S––unday', 'Saturday'), ('S–u–nday', 'Saturday'), ('Su––nday', 'Saturday')]
```

###Hamming distance

This function computes the Hamming distance of two sequences of same length.

Example:

```
>>> from bioinfobox.hamming import hamming
>>> distance = hamming('Slavia', 'Sparta')
>>> print(distance)
3
```

###Fasta processing

The function implemented in bioinfobox.fasta allow for parsing and retreiving
simple information about sequences in the fasta format.

Example:

```
>>> from bioinfobox.fasta import load_fasta, fasta_seq_length
>>> fasta_seqs = load_fasta('python_version/test/fastaseqtest')
>>> record = fasta_seqs['MCHU']
>>> print(fasta_seq_length(record))
150
```

```
>>> from bioinfobox.fasta import load_fasta, fasta_get_seq_and_desc
>>> fasta_seqs = load_fasta('python_version/test/fastaseqtest')
>>> record = fasta_seqs['MCHU']
>>> seq, description = fasta_get_seq_and_desc(record)
>>> print(description)
MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
```

```
>>> from bioinfobox.fasta import load_fasta, fasta_seq_subsequence
>>> fasta_seqs = load_fasta('python_version/test/fastaseqtest')
>>> record = fasta_seqs['MCHU']
>>> subseq = fasta_seq_subsequence(rec, 0, 5)
>>> print(subseq)
MADQL
```

###Multiple sequence alignment

Functions implemented in bioinfobox.msa allow user to parse files 
in the CLUSTAL format, retrieve sequences or columns in the MSA, 
and to compute the sum of pairs score of chosen column or the whole MSA.

Examples:

```
>>> from bioinfobox.msa import load_clustal, clustal_retrieve_position
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> record1 = clustal_retrieve_position(clustal_list, 54)
>>> print(len(record1))
530
```

```
>>> from bioinfobox.msa import *
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> record1 = clustal_retrieve_position(clustal_list, 54)
>>> record2 = clustal_retrieve_id(clustal_list, 'UniRef90_G3WS63')
>>> print(record1 == record2)
True
```

```
>>> from bioinfobox.msa import load_clustal, clustal_column
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> column = clustal_column(clustal_list, 529)
>>> print(''.join(column))
DEEEDDDDDDDDDDDDGDEDDDDDDDDDDDDDDDDDDDDDDNDDDDDDDDDDDDD
```

For scoring a substitution matrix will be necessary. These can be loaded from
bioinfobox.matrix.

Scoring examples:

```
>>> from bioinfobox.msa import load_clustal, sum_of_column
>>> from bioinfobox.matrix import load_matrix
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> column = clustal_column(clustal_list, 529)
>>> mx = load_matrix('BLOSUM62')
>>> print(sum_of_column(column, mx))
7470.0
```

```
>>> from bioinfobox.msa import load_clustal, msa_score
>>> from bioinfobox.matrix import load_matrix
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> mx = load_matrix('BLOSUM62')
>>> print(msa_score(clustal_list, mx))
83756.0
```

###Determinig conservation

The functions in bioinfobox.conservation serve to compute the conservation score
of column or finding the n most conserved positions.

Examples:

```
>>> from bioinfobox.msa import load_clustal
>>> from bioinfobox.conservation import conservation_at_position
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> print(conservation_at_position(clustal_list, 529))
0.8909090909090909
```

```
>>> from bioinfobox.msa import load_clustal
>>> from bioinfobox.conservation import n_most_conserved
>>> clustal_list = load_clustal('python_version/test/p53_mafft_clustal.txt')
>>> print(n_most_conserved(clustal_list, 5))
[(210, 1.0), (211, 1.0), (214, 1.0), (216, 1.0), (219, 1.0)]
```

###Parsing PDB files

The functions in bioinfobox.pdb allow the user to parse a PDB format file and to
retrieve some infomations about the structures.

Examples:

```
>>> from bioinfobox.pdb import load_pdb, pdb_info
>>> pdb_structure = load_pdb('python_version/test/pdb1lb5.ent.txt', '1lb5')
>>> structure_info = pdb_info(pdb_structure)
>>> print(structure_info)
{'models': 1, 'chains': 2, 'residues': 242, 'atoms': 1409}
```

```
>>> from bioinfobox.pdb import load_pdb, pdb_width
>>> pdb_structure = load_pdb('python_version/test/pdb1lb5.ent.txt', '1lb5')
>>> print(pdb_width(pdb_structure))
53.48979
```

Moreover one can obtain the substructure desired with functions pdb_get_models,
pdb_get_chains, pdb_get_residues, and pdb_get_atoms.

Example:

```
>>> from bioinfobox.pdb import load_pdb, pdb_get_chains
>>> pdb_structure = load_pdb('python_version/test/pdb1lb5.ent.txt', '1lb5')
>>> chains = list(pdb_get_chains(pdb_structure))
>>> print(chains)
[<Chain id=A>, <Chain id=B>]
```

Lastly, the user can find neighborhood of a given atom in a structure with
functions pdb_get_neighbor_atoms and pdb_get_neighbor_residues.

Examples:

```
>>> from bioinfobox.pdb import load_pdb, pdb_get_neighbor_atoms
>>> pdb_structure = load_pdb('python_version/test/pdb1lb5.ent.txt', '1lb5')
>>> neigh = pdb_get_neighbor_atoms(pdb_structure[0]['B'][('W', 821, ' ')]['O'], 5, pdb_structure)
>>> print(neigh)
[<Atom O>, <Atom OH>, <Atom CE2>, <Atom CD2>, <Atom CZ>, <Atom CE1>]
```

###Structure and sequence correspondence

Lastly, the function contained in bioinfobox.str_seq allows to analyze the relative
conservation of active site in a structure.

Example:

```
>>> from bioinfobox.str_seq import structure_sequence_equivalence
>>> from bioinfobox.msa import load_clustal
>>> from bioinfobox.pdb import load_pdb
>>> pdb_str = load_pdb('python_version/test/pdb3kz8.ent.txt', '3kz8')
>>> msa = load_clustal('python_version/test/p53_mafft_clustal (kopie).txt')
>>> active_site, average = structure_sequence_equivalence(pdb_str[0]['B'][('H_ZN', 1, ' ')]['ZN'], msa, pdb_str, 80)
>>> print(active_site, average)
0.8363095238095238 0.5388140161725068
```

*In testing this feature I only had access to partial sequence, to which I redsponded by padding the sequence
in the MSA, which may in turn slightly lower the average conservation score, but, I believe, not enough
the compromise the result that shows that the active site is indeed significantly more conserved than
the totality of the structure.*


##Console functionality

This library has also a command line interface, albeit somewhat limited, as I found the biopython
atom format quite impractical to use in this way and so the str_seq and majority of pdb libraries
are not implemented. The rest works as follows:

###Edit distance:

Usage:

```console
$ bioinfobox edit seq1 seq2 [--align]
```

Example:

```console
$ bioinfobox edit Sunday Saturday --align
The edit distance is: 3

S––unday
Saturday

S–u–nday
Saturday

Su––nday
Saturday
```

###Hamming distance:

Usage:

```console
$ bioinfobox hamming seq1 seq2
```

Example:

```console
$ bioinfobox hamming Slavia Sparta
The Hamming distance of the sequences is: 3
```

###Fasta processing:

Usage:

```console
$ bioinfobox fasta path [--sequence id] [--description id] [--len id] [--subsequence id start end]
```

Example:

```console
$ bioinfobox fasta python_version/test/fastaseqtest --description MCHU
MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
```

###Multiple sequence alignment:

Usage:

```console
$ bioinfobox msa path [--position number] [--id id] [--column num] [--score mx_name]
```

Example:

```console
$ bioinfobox msa python_version/test/p53_mafft_clustal.txt --score BLOSUM62
83756.0
```

###Conservation:

Usage:

```console
$ bioinfobox conservation path [--position num] [--best n]
```

Example:

```console
$ bioinfobox conservation python_version/test/p53_mafft_clustal.txt --best 5
[(210, 1.0), (211, 1.0), (214, 1.0), (216, 1.0), (219, 1.0)]
```

###Parsing PDB files:

Due to the low suitability of the biopython Atom format, only the width computation
is implemented.

Usage:

```console
$ bioinfobox pdb path
```

Example:

```console
$ bioinfobox pdb python_version/test/pdb1lb5.ent.txt
The width of the structure is: 53.48979
```
###Matrices:

For the MSA portion it might be useful to know the list of available matrices, which
is what this command provides.

Usage:

```console
$ bioinfobox matrix
```

Example:

```console
$ bioinfobox matrix
1 BENNER22
2 BENNER6
3 BENNER74
4 BLOSUM45
5 BLOSUM50
6 BLOSUM62
7 BLOSUM80
8 BLOSUM90
9 DAYHOFF
10 FENG
11 GENETIC
12 GONNET1992
13 HOXD70
14 JOHNSON
15 JONES
16 LEVIN
17 MCLACHLAN
18 MDM78
19 NUC.4.4
20 PAM250
21 PAM30
22 PAM70
23 RAO
24 RISLER
25 SCHNEIDER
26 STR
27 TRANS
```
