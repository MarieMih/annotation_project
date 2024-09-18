import sys
from Bio import SeqIO
import pandas as pd

def filter_fasta_by_ids(id_file, fasta_file, output_fasta, col):
    """
    Removes records from a FASTA file whose IDs are listed in the second column of a given CSV file.

    Parameters:
    id_file (str): Path to the CSV file containing IDs in the second column.
    fasta_file (str): Path to the input FASTA file.
    output_fasta (str): Path to the output FASTA file to save filtered records.
    """

    df = pd.read_csv(id_file, sep='\t')

    ids_to_remove = set(df.iloc[:, col].astype(str))

    with open(fasta_file, 'r') as fasta_input, open(output_fasta, 'w') as fasta_output:
        for record in SeqIO.parse(fasta_input, 'fasta'):
            if record.id not in ids_to_remove:
                SeqIO.write(record, fasta_output, 'fasta')

filter_fasta_by_ids(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
