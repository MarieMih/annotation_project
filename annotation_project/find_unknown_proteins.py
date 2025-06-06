import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
sys.path.append(os.path.dirname(__file__))
from create_union_protein_fasta_from_gffs import get_from_upimapi
import common_variables


def get_length_of_proteins(faa_file):
    protein_data = []
    for record in SeqIO.parse(faa_file, "fasta"):
        protein_id = record.id
        protein_length = len(record.seq)
        protein_data.append([protein_id, protein_length])

    df = pd.DataFrame(protein_data, columns=["Locus Tag", "Length"])
    return df


def finding_from_fasta(fasta_file):
    """
    Take faa-file with representative unknown protein for each cluster from all samples
    (from directory "union_unknown_faa" and file "union_rep_seq.fasta")
    and search them in Uniprot with UPIMAPI.
    """
    upimapi_output = os.path.join(os.path.split(fasta_file)[0], "upimapi_unknown_protein_results")
    subprocess.run(["upimapi",
                    '-i', fasta_file,
                    '-o', upimapi_output,
                    '-t', common_variables.N_THREADS,
                    '-rd', common_variables.UPIMAPI_RESOURCES,
                    '-db', common_variables.UPIMAPI_DATABASE
                    ],
                   check=True)

    results = pd.read_csv(os.path.join(upimapi_output, "UPIMAPI_results.tsv"), sep="\t", header=0)

    upimapi_fasta = os.path.join(os.path.split(fasta_file)[0], "upimapi_fasta")

    get_from_upimapi(set(results["Entry"]), upimapi_fasta)
    query_len_df = get_length_of_proteins(fasta_file)

    results['Reference a_length'] = 1
    results['Query a_length'] = -1


    ### filtering results by percent of identity, length of query and length of target

    with open(os.path.join(upimapi_fasta, "upimapi_output", "uniprotinfo.fasta"), 'r') as main_faa:
        for record in SeqIO.parse(main_faa, "fasta"):
            loctag = results.loc[results["Entry"] == record.id.split("|")[1], 'qseqid']
            tmp_len = query_len_df.loc[query_len_df["Locus Tag"] == loctag.values[0], 'Length']

            results.loc[results["Entry"] == record.id.split("|")[1], 'Reference a_length'] = len(record.seq)
            results.loc[results["Entry"] == record.id.split("|")[1], 'Query a_length'] = tmp_len.values[0]

    results['Percent of length'] = results['Reference a_length'] / results['Query a_length']
    results['Percent of align'] = results['Query a_length'] / results['length']
    results = results[(results["pident"] > common_variables.PIDENT)
                      & (results['Percent of length'] > common_variables.RATIO_LEN_BOTTOM)
                      & (results['Percent of length'] < common_variables.RATIO_LEN_UP)
                      & (results['Percent of align'] > common_variables.QUERY_ALIGN_LENGTH_BOTTOM)
                      & (results['Percent of align'] < common_variables.QUERY_ALIGN_LENGTH_UP)]
    return results
