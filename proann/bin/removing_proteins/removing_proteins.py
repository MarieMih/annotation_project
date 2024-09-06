from Bio import SeqIO
import sys

# Read IDs to be removed
with open(sys.argv[3]) as f:
    ids_to_remove = set(line.strip() for line in f)

# Function to check if any of the IDs is in the header
def should_remove(record):
    header = record.description
    return any(id_ in header for id_ in ids_to_remove)

# Read the input FASTA file, filter, and write to the output FASTA file
with open(sys.argv[1]) as input_handle, open(sys.argv[2], "w") as output_handle:
    records = (record for record in SeqIO.parse(input_handle, "fasta") if not should_remove(record))
    SeqIO.write(records, output_handle, "fasta")
