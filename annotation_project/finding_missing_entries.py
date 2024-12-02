import os
import subprocess
import pandas as pd
from Bio import SeqIO
from create_union_protein_fasta_from_gffs import get_from_upimapi

UPIMAPI_RESOURCES = "/storage/data1/marmi/upimapi_databases"
UPIMAPI_DATABASE = "uniprot"  # change on uniprot

QUERY_ALIGN_LENGTH_UP = 1.05
QUERY_ALIGN_LENGTH_BOTTOM = 0.95
PIDENT = 90
RATIO_LEN_UP = 1.1
RATIO_LEN_BOTTOM = 0.9


def finding_missing_entries(df, faa):
    """
    Take all record with NaN in Entry UniProtKB column and annotate it with upimapi.
    """

    missing_entries = df[df["Entry UniProtKB"].isna()]
    missing_entries = missing_entries[missing_entries["Type"] == "cds"]
    missing_entries['A_length'] = abs(missing_entries["Stop"] - missing_entries["Start"] + 1) // 3
    missing_entries_set = set(missing_entries["Locus Tag"])

    new_faa_file = faa.rpartition(".")[0] + "_missing_KB.faa"

    with open(new_faa_file, 'w') as recordfile:
        with open(faa, 'r') as main_faa:
            for record in SeqIO.parse(main_faa, "fasta"):
                if record.id in missing_entries_set:
                    SeqIO.write(record, recordfile, "fasta")

    upimapi_output = os.path.join(os.path.split(new_faa_file)[0], "upimapi_missing_results")
    subprocess.run(["upimapi",
                    '-i', new_faa_file,
                    '-o', upimapi_output,
                    '-t', '8',
                    '-rd', UPIMAPI_RESOURCES,
                    '-db', UPIMAPI_DATABASE
                    ],
                   check=True)

    results = pd.read_csv(os.path.join(upimapi_output, "UPIMAPI_results.tsv"), sep="\t", header=0)

    upimapi_fasta = os.path.join(os.path.split(new_faa_file)[0], "upimapi_fasta")

    get_from_upimapi(set(results["Entry"]), upimapi_fasta)

    results['Reference a_length'] = 1
    results['Query a_length'] = -1

    with open(os.path.join(upimapi_fasta, "upimapi_output", "uniprotinfo.fasta"), 'r') as main_faa:
        for record in SeqIO.parse(main_faa, "fasta"):
            loctag = results.loc[results["Entry"] == record.id.split("|")[1], 'qseqid']
            tmp_len = missing_entries.loc[missing_entries["Locus Tag"] == loctag.values[0], 'A_length']

            results.loc[results["Entry"] == record.id.split("|")[1], 'Reference a_length'] = len(record.seq)
            results.loc[results["Entry"] == record.id.split("|")[1], 'Query a_length'] = tmp_len.values[0]

    results['Percent of length'] = results['Reference a_length'] / results['Query a_length']
    results['Percent of align'] = results['Query a_length'] / results['length']

    results = results[(results["pident"] > PIDENT)
                      & (results['Percent of length'] > RATIO_LEN_BOTTOM)
                      & (results['Percent of length'] < RATIO_LEN_UP)
                      & (results['Percent of align'] > QUERY_ALIGN_LENGTH_BOTTOM)
                      & (results['Percent of align'] < QUERY_ALIGN_LENGTH_UP)]

    for i, r in results.iterrows():

        # information from UPIMAPI_results.tsv
        ltag = r["qseqid"]
        # organism = r["Taxonomic lineage (SPECIES)"]
        # gene = str(r["Gene Names"]).split(" ")[0]
        # entry = r["sseqid"]
        # product = r["Protein names"]

        # # write this information to original df
        # df.loc[df['Locus Tag'] == ltag, ["Organism"]] = organism
        # df.loc[df['Locus Tag'] == ltag, ["Gene"]] = gene
        # df.loc[df['Locus Tag'] == ltag, ["Entry UniProtKB"]] = entry
        # df.loc[df['Locus Tag'] == ltag, ["Product"]] = product

        # write this information to original df
        df.loc[df['Locus Tag'] == ltag, ["Organism"]] = r["Taxonomic lineage (SPECIES)"]
        df.loc[df['Locus Tag'] == ltag, ["Gene"]] = str(r["Gene Names"]).split(" ")[0]
        df.loc[df['Locus Tag'] == ltag, ["Entry UniProtKB"]] = r["sseqid"]
        df.loc[df['Locus Tag'] == ltag, ["Product"]] = r["Protein names"]
        df.loc[df['Locus Tag'] == ltag, ['GO']] = r['Gene Ontology (GO)']  # new!!!!!!!!!!
        df.loc[df['Locus Tag'] == ltag, ['KEGG']] = r['KEGG']  # new!!!!!!!!!!
        df.loc[df['Locus Tag'] == ltag, ['UniPathway']] = r['UniPathway']  # new!!!!!!!!!!
        df.loc[df['Locus Tag'] == ltag, ['Pathway']] = r['Pathway']  # new!!!!!!!!!!
        df.loc[df['Locus Tag'] == ltag, ['Keywords']] = r['Keywords']  # new!!!!!!!!!!
