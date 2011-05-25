#!/usr/bin/python

from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment

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

for id, clade in enumerate(tree.find_clades()):
    clade.our_unique_identifier = id

# this isn't quite done yet.....
# we need to figure out the issue of assigning unique ids 
# to internal nodes
def build_msa(node, sequence_msa_map):
    key = str(id(node))
    file = key + '.msa'
    file_handle = open(file, 'w')
    terminals = node.get_terminals()
    alignments = [sequence_msa_map[terminal.name] for terminal in
            terminals]
    alignments = MultipleSeqAlignment(alignments)
    AlignIO.write(alignments, file_handle, 'stockholm')

build_msa(tree.root, sequence_msa_map)
