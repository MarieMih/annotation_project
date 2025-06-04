import csv
from Bio import SeqIO


def divide_fasta(input_file):
    """
    This function gets input annotation file from bakta (.tsv)
    and divides into six files.
    """

    fasta_file = input_file.replace(".tsv", ".faa")

    user_protein_set = set()
    gene_symbol_set = set()
    unknown_set = set()

    prefix = input_file.rpartition('.')[0]
    user_protein_fasta = prefix + "_userproteins.faa"
    gene_detected_fasta = prefix + "_detected.faa"
    unknown_fasta = prefix + "_unknown.faa"

    with open(input_file.replace(".tsv", "_extended.tsv"), 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if row[0].startswith("#"):
                continue
            if (row[1] != "cds") and (row[1] != "sorf"):  # column "Type"
                continue
            row_str = '\t'.join(row)
            if 'UserProtein' in row_str:
                user_protein_set.add(row[5])  # column "Locus Tag"
            elif ('UserProtein' not in row_str) and (row[6] != ""):  # column "Gene"
                gene_symbol_set.add(row[5])
            else:
                unknown_set.add(row[5])

    user_protein_records = list()
    gene_detected_records = list()
    unknown_records = list()

    fasta_records = SeqIO.parse(fasta_file, "fasta")
    for record in fasta_records:
        if record.id in user_protein_set:
            user_protein_records.append(record)
        elif record.id in gene_symbol_set:
            gene_detected_records.append(record)
        elif record.id in unknown_set:
            unknown_records.append(record)

    with open(user_protein_fasta, 'w') as output_handle:
        SeqIO.write(user_protein_records, output_handle, "fasta")
    with open(gene_detected_fasta, 'w') as output_handle:
        SeqIO.write(gene_detected_records, output_handle, "fasta")
    with open(unknown_fasta, 'w') as output_handle:
        SeqIO.write(unknown_records, output_handle, "fasta")

    print("Protein fasta parsed.")
