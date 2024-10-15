import os
import subprocess
import pandas as pd
from Bio import SeqIO

UPIMAPI_RESOURCES = "/storage/data1/marmi/upimapi_databases"
UPIMAPI_DATABASE = "swissprot"


def finding_missing_entries(df, faa):
    """
    Take all record with NaN in Entry UniProtKB column and annotate it with upimapi.
    """

    missing_entries = df[df["Entry UniProtKB"].isna()]
    missing_entries = missing_entries[missing_entries["Type"] == "cds"]
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
    results = results[results["pident"] > 90]

    for i, r in results.iterrows():

        # information from UPIMAPI_results.tsv
        ltag = r["qseqid"]
        organism = r["Taxonomic lineage (SPECIES)"]
        gene = r["Gene Names"].split(" ")[0]
        entry = r["sseqid"]
        product = r["Protein names"]

        # write this information to original df
        df.loc[df['Locus Tag'] == ltag, ["Organism"]] = organism
        df.loc[df['Locus Tag'] == ltag, ["Gene"]] = gene
        df.loc[df['Locus Tag'] == ltag, ["Entry UniProtKB"]] = entry
        df.loc[df['Locus Tag'] == ltag, ["Product"]] = product

