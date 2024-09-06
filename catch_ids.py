import csv
from Bio import SeqIO
import sys

def extract_fasta_records(tsv_file, fasta_file, output_fasta_file):
    # Read the list of IDs from the TSV file (6th column)
    ids = set()
    with open(tsv_file, 'r') as tsv_handle:
        reader = csv.reader(tsv_handle, delimiter='\t')
        for row in reader:
            if len(row) > 5:  # Ensure there is a 6th column
                ids.add(row[5])
    
    # Parse the FASTA file and extract the records with matching IDs
    fasta_records = SeqIO.parse(fasta_file, "fasta")
    matching_records = [record for record in fasta_records if record.id in ids]

    # Write the matching records to the output FASTA file
    with open(output_fasta_file, 'w') as output_handle:
        SeqIO.write(matching_records, output_handle, "fasta")

    print(f"Extracted {len(matching_records)} records to {output_fasta_file}")

# Example usage
tsv_file = sys.argv[1]  # Path to your TSV file
fasta_file = sys.argv[2]  # Path to your FASTA file
output_fasta_file = tsv_file.rpartition('.')[0]+"_by_bakta_tag.faa"  # Path to the output FASTA file

extract_fasta_records(tsv_file, fasta_file, output_fasta_file)
