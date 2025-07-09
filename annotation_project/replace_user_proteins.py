import subprocess
import os
import re
import csv
import pandas as pd


def get_user_protein_information(tsv):
    """
    Get information about UserProtein id.
    Input: '_userproteins_only.tsv'
    Output: directory from upimapi
    """
    pth = os.path.split(os.path.abspath(tsv))[0]

    list_ids = os.path.join(pth, "userproteins_ids.csv")
    upimapi_dir = os.path.join(pth, "userprotein_upimapi")

    uniprotkb_ids = set()

    df_origin = pd.DataFrame(columns=["Locus Tag", "UserProtein"])

    with open(tsv, 'r') as f:

        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            row_str = '\t'.join(row)
            uniprotkb = row_str
            entry = re.compile(r"UserProtein:UniProtKB\|(\w+)")
            if isinstance(uniprotkb, str):
                match = re.search(entry, uniprotkb)
                if match:
                    uniprotkb = match.group(1)
                    uniprotkb_ids.add(uniprotkb)
                    df_origin.loc[len(df_origin)] = {'Locus Tag': row[5], 'UserProtein': uniprotkb}

    with open(list_ids, 'w') as f:
        f.write(",".join(uniprotkb_ids))

    subprocess.run(['upimapi', '-i', list_ids,
                    '-o', upimapi_dir,
                    '--columns', "Entry&Entry Name&Gene Names&Protein names&EC number&Function [CC]&Pathway&Keywords&Protein existence&Gene Ontology (GO)&Protein families&Taxonomic lineage&Taxonomic lineage (Ids)&Taxonomic lineage IDs (SPECIES)&Taxonomic lineage (SPECIES)&Organism&Organism (ID)&BioCyc&BRENDA&CDD&eggNOG&Ensembl&InterPro&KEGG&Pfam&Reactome&RefSeq&UniPathway",
                    '-t', '1'],
                   check=True)

    df = pd.read_csv(os.path.join(upimapi_dir, "uniprotinfo.tsv"), sep="\t", header=0)
    df['Gene Names'] = df['Gene Names'].str.split().str[0]

    df = df[["Entry", "Gene Names", "Protein names", "Organism", "Gene Ontology (GO)", "KEGG", "UniPathway", "Pathway", "Keywords"]]
    df = df.merge(df_origin, how='right', left_on='Entry', right_on='UserProtein')
    df = df.drop(columns=['UserProtein'])

    return df
