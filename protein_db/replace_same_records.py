import sys
from Bio import SeqIO


def replace_unique_fasta(input_file, output_file):
    all_records = []
    unique_headers = []
    header_keys = []

    # Parse the FASTA file and replace unique headers
    for record in SeqIO.parse(input_file, "fasta"):
        header_key = " ".join(record.description.split(" ")[1:])
        if header_key not in header_keys:
            all_records.append(record)
            unique_headers.append(record.description)
            header_keys.append(header_key)
        else:
            new_header = [i for i in unique_headers if header_key in i]
            assert new_header[0] != ''
            record.description = new_header[0]
            record.id = new_header[0].split(' ')[0]
            all_records.append(record)

    with open(output_file, "w") as output_handle:
        SeqIO.write(all_records, output_handle, "fasta")


input_file = sys.argv[1]
output_file = sys.argv[2]
replace_unique_fasta(input_file, output_file)
