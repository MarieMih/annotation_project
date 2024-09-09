from Bio import SeqIO
import sys

def filter_unique_fasta(input_file, output_file):
    unique_records = {}
    unique_headers = []

    # Parse the FASTA file and filter unique headers
    for record in SeqIO.parse(input_file, "fasta"):
        # Split only once on space and check uniqueness using the first part of the header
        header_key = " ".join(record.description.split(" ")[1:])
        if header_key not in unique_headers:
            unique_records[record.id] = record
            unique_headers.append(header_key)

    # Write the unique records to the output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(unique_records.values(), output_handle, "fasta")

# Usage example
input_file = sys.argv[1]
output_file = sys.argv[2]
filter_unique_fasta(input_file, output_file)
