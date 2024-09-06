import os, sys
from Bio import SeqIO

my_file = sys.argv[1]  # Obviously not FASTA
print(my_file)
def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

print(is_fasta(my_file))

