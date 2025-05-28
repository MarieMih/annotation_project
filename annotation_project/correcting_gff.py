import os
import sys
import shutil
import re
import csv
import pandas as pd
sys.path.append(os.path.dirname(__file__))
from replace_user_proteins import get_user_protein_information


def correcting_gff(input_path):
    """
    creating gff and tsv with modified information
    """

    ###################### declaration of variables ######################

    ######################################################################
    #  All analysis divides into two parts: unknown proteins and         #
    #  proteins with known uniref100 references.                         #
    ######################################################################

    bakta_gff = ""

    for file in os.listdir(input_path):
        if file.endswith(".gff3"):
            bakta_gff = os.path.abspath(input_path + "/" + file)
            bakta_tsv = os.path.abspath(input_path + "/" + file.rpartition('.')[0] + ".tsv")

    if os.path.exists(bakta_gff):
        print(f'The file {bakta_gff} exists')
    else:
        print(f'The file {bakta_gff} does not exist')
        sys.exit()

    bakta_gff_ext = bakta_gff.rpartition('.')[0] + "_extended.gff3"
    bakta_tsv_ext = bakta_tsv.rpartition('.')[0] + "_extended.tsv"

    faa = bakta_gff.rpartition('.')[0] + ".faa"

    shutil.copy(bakta_gff, bakta_gff_ext)
    shutil.copy(bakta_tsv, bakta_tsv_ext)

    pth = os.path.split(os.path.abspath(bakta_gff_ext))[0]

    for file in os.listdir(pth):
        if file.endswith("ref2ref"):
            uni = os.path.join(os.path.abspath(pth), file)  # find directory with output of upimapi

    for file in os.listdir(pth):
        if file.endswith("_uniref100_columns.tsv"):
            # in this file first column is a unique protein id from bakta and a second column is a uniref100 reference from bakta
            uniref100_table_locus_tag_tsv = os.path.abspath(pth) + "/" + file
        if file.endswith("_userproteins_only.tsv"):
            userproteins_only = os.path.join(os.path.abspath(pth), file)

    # this file contains all information from uniprot for found representative entry of known part bakta annotation
    uniref100_representative_tsv = uni + "/uniprotkb/uniprotinfo.tsv"

    ############################### meaningful part ###########################

    extended_tsv_df = pd.read_csv(bakta_tsv_ext, sep="\t",
                                  header=None,
                                  comment='#',
                                  names=['Sequence Id', 'Type', 'Start', 'Stop', 'Strand', 'Locus Tag', 'Gene', 'Product', 'DbXrefs'])
    extended_tsv_df["Organism"] = "Escherichia coli"
    extended_tsv_df.to_csv(bakta_tsv_ext, sep='\t', index=False)

    ### known proteins

    joined_uniref_df = pd.read_csv(uniref100_table_locus_tag_tsv, sep="\t", header=None)
    joined_uniref_df.columns = ['id', 'uniprot']

    uniref_df_part = pd.read_csv(uniref100_representative_tsv, sep="\t", header=0)
    joined_uniref_df = joined_uniref_df.merge(uniref_df_part, left_on="uniprot", right_on="Entry", how='left', suffixes=('', '_new'))
    joined_uniref_df['Gene Names'] = joined_uniref_df['Gene Names'].str.split().str[0]

    extended_tsv_df = pd.read_csv(bakta_tsv_ext, sep="\t", header=0)
    merged_df = extended_tsv_df.merge(joined_uniref_df, left_on="Locus Tag", right_on='id', suffixes=('', '_new'), how="left")
    joined_uniref_df.fillna('', inplace=True)


####### begin - rewrite UniRef100 records
    for i in range(len(merged_df)):

        locus_tag = merged_df.at[i, 'Locus Tag']
        matching_row = joined_uniref_df[joined_uniref_df['id'] == locus_tag]
        if not matching_row.empty:
            if str(matching_row['Gene Names'].values[0]) != "":
                merged_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
            if str(matching_row['Protein names'].values[0]) != "":
                merged_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
            merged_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
            merged_df.at[i, 'Entry UniProtKB'] = matching_row['Entry'].values[0]
            merged_df.at[i, 'GO'] = matching_row['Gene Ontology (GO)'].values[0]
            merged_df.at[i, 'KEGG'] = matching_row['KEGG'].values[0]
            merged_df.at[i, 'UniPathway'] = matching_row['UniPathway'].values[0]
            merged_df.at[i, 'Pathway'] = matching_row['Pathway'].values[0]
            merged_df.at[i, 'Keywords'] = matching_row['Keywords'].values[0]

    extended_tsv_df["Gene"] = merged_df["Gene"]
    extended_tsv_df["Product"] = merged_df["Product"]
    extended_tsv_df["Organism"] = merged_df["Organism"]
    extended_tsv_df["Entry UniProtKB"] = merged_df["Entry UniProtKB"]
    extended_tsv_df["GO"] = merged_df["GO"]
    extended_tsv_df["KEGG"] = merged_df["KEGG"]
    extended_tsv_df['UniPathway'] = merged_df['UniPathway']
    extended_tsv_df['Pathway'] = merged_df['Pathway']
    extended_tsv_df['Keywords'] = merged_df['Keywords']
####### end - rewrite UniRef100 records

####### begin - rewrite fields with UserProtein data
    try:
        userprotein_df = get_user_protein_information(userproteins_only)

        for i in range(len(extended_tsv_df)):

            locus_tag = extended_tsv_df.at[i, 'Locus Tag']
            matching_row = userprotein_df[userprotein_df['Locus Tag'] == locus_tag]

            if not matching_row.empty:

                if str(matching_row['Gene Names'].values[0]) != "":
                    extended_tsv_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
                if str(matching_row['Protein names'].values[0]) != "":
                    extended_tsv_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
                extended_tsv_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
                extended_tsv_df.at[i, 'Entry UniProtKB'] = matching_row['Entry'].values[0]
                extended_tsv_df.at[i, 'GO'] = matching_row['Gene Ontology (GO)'].values[0]
                extended_tsv_df.at[i, 'KEGG'] = matching_row['KEGG'].values[0]
                extended_tsv_df.at[i, 'UniPathway'] = matching_row['UniPathway'].values[0]
                extended_tsv_df.at[i, 'Pathway'] = matching_row['Pathway'].values[0]
                extended_tsv_df.at[i, 'Keywords'] = matching_row['Keywords'].values[0]
    except:
        print("Userprotein_only file is empty.")
####### end - rewrite fields with UserProtein data

####### begin - fill empty fields with "", add Transcript_id/Gene_id
    extended_tsv_df["Gene"] = extended_tsv_df["Gene"].fillna('').astype(str)
    extended_tsv_df["Entry UniProtKB"] = extended_tsv_df["Entry UniProtKB"].fillna('').astype(str)
    extended_tsv_df["Transcript_id"] = extended_tsv_df["Type"] + "|" + extended_tsv_df["Gene"].fillna('') + "|" + extended_tsv_df["Entry UniProtKB"].fillna('')
    for i in range(len(extended_tsv_df)):
        if ((extended_tsv_df.at[i, 'Gene'] == "")
           and (extended_tsv_df.at[i, 'Entry UniProtKB'] == "")
           and ((extended_tsv_df.at[i, 'Type'] == "cds") or (extended_tsv_df.at[i, 'Type'] == "sorf"))):
            extended_tsv_df.at[i, 'Transcript_id'] = str(extended_tsv_df.at[i, "Locus Tag"]) + "_" + str(extended_tsv_df.at[i, "Type"])
    extended_tsv_df["Gene_id"] = extended_tsv_df["Transcript_id"]
####### end - fill empty fields with "", add Transcript_id/Gene_id

####### begin - save to tsv and to gff
    extended_tsv_df.to_csv(bakta_tsv_ext, sep='\t', index=False)
    extended_tsv_df.fillna('', inplace=True)

    ### creating new gff

    with open(bakta_gff_ext, 'w') as w:
        with open(bakta_gff, 'r') as f:

            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row[0].startswith('#'):
                    w.write('\t'.join(row) + "\n")
                    continue
                if len(row) < 2:
                    w.write('\t'.join(row) + "\n")
                    continue
                pairs = row[8].split(';')
                parsed_dict = {pair.split('=', 1)[0]: pair.split('=', 1)[1] for pair in pairs}
                id = extended_tsv_df[extended_tsv_df['Locus Tag'] == parsed_dict.get('locus_tag')]
                if id.empty:
                    w.write('\t'.join(row)+"\n")
                else:
                    record = id.iloc[0]
                    new_record = []
                    if record['Locus Tag']:
                        new_record.append(f"ID={record['Locus Tag']}")
                    if record['Product']:
                        new_record.append(f"Name={record['Product'].replace(';', ',')}")
                    if record['Locus Tag']:
                        new_record.append(f"locus_tag={record['Locus Tag']}")
                    if record['Product']:
                        new_record.append(f"product={record['Product'].replace(';', ',')}")
                    if record['DbXrefs']:
                        new_record.append(f"Dbxref={record['DbXrefs']}")
                    if record['Gene']:
                        new_record.append(f"gene={record['Gene']}")
                    if record['Entry UniProtKB']:
                        new_record.append(f"entry={record['Entry UniProtKB']}")
                    if record['Organism']:
                        new_record.append(f"organism={record['Organism']}")
                    if record['Transcript_id']:
                        new_record.append(f"transcript_id={record['Transcript_id']}")
                    if record['Gene_id']:
                        new_record.append(f"gene_id={record['Gene_id']}")
                    if record['GO']:
                        new_record.append(f"go={record['GO']}")
                    if record['KEGG']:
                        new_record.append(f"kegg={record['KEGG']}")
                    if record['UniPathway']:
                        new_record.append(f"unipathway={record['UniPathway']}")
                    if record['Pathway']:
                        new_record.append(f"pathway={record['Pathway']}")
                    if record['Keywords']:
                        new_record.append(f"keywords={record['Keywords']}")

                    new_record = ";".join(new_record)
                    new_row = row.copy()
                    new_row[8] = new_record
                    w.write('\t'.join(new_row)+"\n")
####### end - save to tsv and to gff

    return bakta_gff_ext
