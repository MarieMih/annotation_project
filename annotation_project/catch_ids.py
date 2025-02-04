import csv
from Bio import SeqIO


def extract_fasta_records(tsv_file, fasta_file, output_fasta_file):
    """
    Read the list of IDs from the TSV file (6th column)
    """

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


def catch_ids(tsv_file, fasta_file):
    """
    Take '_cds_sorf.tsv' (that contains locus_tag of records without bakta-UniRef100 or UserProteins-UniProtKB)
    and find all UniProt IDs.
    """
    output_fasta_file = tsv_file.rpartition('.')[0]+"_by_bakta_tag.faa"
    extract_fasta_records(tsv_file, fasta_file, output_fasta_file)
