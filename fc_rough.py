#!/usr/bin/python

from Bio import Phylo, AlignIO

tree_file = "sample_tree.ml"
msa_file = "sample_alignment.msa"
new_sequence_file = "new_sequence.fasta"

tree = Phylo.read(tree_file, "newick")
msa = AlignIO.read(msa_file, "fasta")

# Note: the name of the tree leaf nodes is going to be the 
# same as the id of each msa entry

sequence_msa_map = {}
for entry in msa:
    sequence_msa_map[entry.id] = entry


