import sys
from Bio import SeqIO


def replace_unique_fasta(input_file, output_file):
    unique_headers = {}

    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            header_key = " ".join(record.description.split(" ")[1:])

            if header_key not in unique_headers:
                unique_headers[header_key] = (record.description, record.id)
                SeqIO.write(record, output_handle, "fasta")
            else:
                existing_description, existing_id = unique_headers[header_key]
                record.description = existing_description
                record.id = existing_id
                SeqIO.write(record, output_handle, "fasta")


input_file = sys.argv[1]
output_file = sys.argv[2]
replace_unique_fasta(input_file, output_file)
