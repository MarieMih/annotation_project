import sys
from Bio import SeqIO


def filter_unique_fasta(in_file, out_file):
    unique_headers = set()

    with open(out_file, "w") as output_handle:
        for record in SeqIO.parse(in_file, "fasta"):
            header_key = " ".join(record.description.split(" ")[1:])
            if header_key not in unique_headers:
                unique_headers.add(header_key)
                SeqIO.write(record, output_handle, "fasta")


input_file = sys.argv[1]
output_file = sys.argv[2]
filter_unique_fasta(input_file, output_file)
