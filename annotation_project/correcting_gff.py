import os
import shutil
import re
import csv
import pandas as pd
from finding_missing_entries import finding_missing_entries
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

    for file in os.listdir(input_path):
        if file.endswith(".gff3"):
            bakta_gff = os.path.abspath(input_path + "/" + file)
            bakta_tsv = os.path.abspath(input_path + "/" + file.rpartition('.')[0] + ".tsv")

    if os.path.exists(bakta_gff):
        print(f'The file {bakta_gff} exists')
    else:
        print(f'The file {bakta_gff} does not exist')
        exit()

    bakta_gff_ext = bakta_gff.rpartition('.')[0] + "_extended.gff3"
    bakta_tsv_ext = bakta_tsv.rpartition('.')[0] + "_extended.tsv"

    faa = bakta_gff.rpartition('.')[0] + ".faa"

    shutil.copy(bakta_gff, bakta_gff_ext)
    shutil.copy(bakta_tsv, bakta_tsv_ext)

    pth = os.path.split(os.path.abspath(bakta_gff_ext))[0]

    for file in os.listdir(pth):
        if file.endswith("ref2ref"):
            uni = os.path.abspath(pth) + "/" + file  # find directory with output of upimapi

    for file in os.listdir(pth):
        if file.endswith("_uniref100_columns.tsv"):
            # in this file first column is a unique protein id from bakta and a second column is a uniref100 reference from bakta
            info_uniref100_table_tsv = os.path.abspath(pth) + "/" + file
            # in this file first column is a uniref100 reference and a second column is a representative member id entry of cluster
            info_uniref100_table_ids = uni + "/uniprotinfo_uniref_representative_ids.tsv"
        if file.endswith("_userproteins_only.tsv"):
            userproteins_only = os.path.join(os.path.abspath(pth), file)

    # this file contains all information from uniprot for found representative entry of known part bakta annotation
    info_uniref100_tsv = uni + "/uniprotkb/uniprotinfo.tsv"
    # this file contains all information from uniprot for found representative entry of unknown part bakta annotation
    info_unknown_tsv = uni + "/uniprotkb/annotation/uniprotinfo.tsv"
    # in this file first column is a unique protein id from bakta and a second column is a representative member id entry of cluster (unknown part)
    info_unknown_table_tsv = uni + "/uniprotkb/annotation/UPIMAPI_results.tsv"

    unknown_df = pd.read_csv(info_unknown_table_tsv, sep="\t", usecols=range(2), header=0)
    unknown_info = pd.read_csv(info_unknown_tsv, sep="\t", header=0, index_col=0)

    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', 30)

    ############################### meaningful part ###########################

    ext_tsv_df = pd.read_csv(bakta_tsv_ext, sep="\t", header=None, comment='#',
                             names=['Sequence Id', 'Type', 'Start', 'Stop', 'Strand', 'Locus Tag', 'Gene', 'Product', 'DbXrefs'])

    ### unknown proteins

    joined_unknown = unknown_df.join(unknown_info, on="sseqid")
    joined_unknown['Gene Names'] = joined_unknown['Gene Names'].str.split().str[0]  # could be a lot of gene names, taking only first
    merged_df = ext_tsv_df.merge(joined_unknown, left_on="Locus Tag", right_on='qseqid', suffixes=('', '_new'), how="left")
    joined_unknown.fillna('', inplace=True)

    for i in range(len(merged_df)):

        locus_tag = merged_df.at[i, 'Locus Tag']
        matching_row = joined_unknown[joined_unknown['qseqid'] == locus_tag]

        if not matching_row.empty and matching_row['Gene Names'].values[0] != '':
            merged_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
            merged_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
            merged_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
            merged_df.at[i, 'Entry UniProtKB'] = matching_row['sseqid'].values[0]

    ext_tsv_df["Gene"] = merged_df["Gene"]
    ext_tsv_df["Product"] = merged_df["Product"]
    ext_tsv_df["Organism"] = merged_df["Organism"]
    ext_tsv_df.to_csv(bakta_tsv_ext, sep='\t', index=False)

    ### known proteins

    uniref_df_1 = pd.read_csv(info_uniref100_table_tsv, sep="\t", header=None)
    uniref_df_2 = pd.read_csv(info_uniref100_table_ids, sep="\t", header=None)
    uniref_df_2 = uniref_df_2.drop_duplicates(subset=0, keep='first')

    joined_uni = uniref_df_1.merge(uniref_df_2, left_on=1, right_on=0, how='left', suffixes=('', '_new'), validate='many_to_one')
    joined_uni = joined_uni[["0", "1_new"]]
    joined_uni.columns = ['id', 'uniprot']

    uniref_df_3 = pd.read_csv(info_uniref100_tsv, sep="\t", header=0)
    joined_uni = joined_uni.merge(uniref_df_3, left_on="uniprot", right_on="Entry", how='left', suffixes=('', '_new'))
    joined_uni['Gene Names'] = joined_uni['Gene Names'].str.split().str[0]
    ext_tsv_df = pd.read_csv(bakta_tsv_ext, sep="\t", header=0)
    merged_df = ext_tsv_df.merge(joined_uni, left_on="Locus Tag", right_on='id', suffixes=('', '_new'), how="left")
    joined_uni.fillna('', inplace=True)

    for i in range(len(merged_df)):

        locus_tag = merged_df.at[i, 'Locus Tag']
        uniprotkb = merged_df.at[i, 'DbXrefs']
        entry = re.compile(r"UserProtein:[^|]*\|([^,\n]*)")

        matching_row = joined_uni[joined_uni['id'] == locus_tag]

        if not matching_row.empty:

            merged_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
            merged_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
            merged_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
            merged_df.at[i, 'Entry UniProtKB'] = matching_row['Entry'].values[0]

        # if isinstance(uniprotkb, str):  # v2
        #     match = re.search(entry, uniprotkb)
        #     if match:
        #         uniprotkb = match.group(1)
        #         merged_df.at[i, 'Entry UniProtKB'] = uniprotkb
        #         merged_df.at[i, 'Organism'] = "Escherichia coli"

    ext_tsv_df["Gene"] = merged_df["Gene"]
    ext_tsv_df["Product"] = merged_df["Product"]
    ext_tsv_df["Organism"] = merged_df["Organism"]
    ext_tsv_df["Entry UniProtKB"] = merged_df["Entry UniProtKB"]



    ### new block

    try:
        userprotein_df = get_user_protein_information(userproteins_only)

        for i in range(len(ext_tsv_df)):

            locus_tag = ext_tsv_df.at[i, 'Locus Tag']

            matching_row = userprotein_df[userprotein_df['Locus Tag'] == locus_tag]

            if not matching_row.empty:

                ext_tsv_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
                ext_tsv_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
                ext_tsv_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
                ext_tsv_df.at[i, 'Entry UniProtKB'] = matching_row['Entry'].values[0]

    except:
        print("Userprotein_only file is empty.")

    ### new block




    finding_missing_entries(ext_tsv_df, faa)

    ext_tsv_df["Transcript_id"] = ext_tsv_df["Type"] + "|" + ext_tsv_df["Gene"].fillna('') + "|" + ext_tsv_df["Entry UniProtKB"].fillna('')
    ext_tsv_df["Gene_id"] = ext_tsv_df["Transcript_id"]

    ext_tsv_df.to_csv(bakta_tsv_ext, sep='\t', index=False)
    ext_tsv_df.fillna('', inplace=True)

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
                id = ext_tsv_df[ ext_tsv_df['Locus Tag'] == parsed_dict.get('locus_tag')]
                if id.empty:
                    w.write('\t'.join(row)+"\n")
                else:
                    record = id.iloc[0]
                    new_record = []
                    if record['Locus Tag']:
                        ID = f"ID={record['Locus Tag']}"
                        new_record.append(ID)
                    if record['Product']:
                        Name = f"Name={record['Product'].replace(';', ',')}"
                        new_record.append(Name)
                    if record['Locus Tag']:
                        locus_tag = f"locus_tag={record['Locus Tag']}"
                        new_record.append(locus_tag)
                    if record['Product']:
                        product = f"product={record['Product'].replace(';', ',')}"
                        new_record.append(product)
                    if record['DbXrefs']:
                        Dbxref = f"Dbxref={record['DbXrefs']}"
                        new_record.append(Dbxref)
                    if record['Gene']:
                        gene = f"gene={record['Gene']}"
                        new_record.append(gene)
                    if record['Entry UniProtKB']:
                        entry = f"entry={record['Entry UniProtKB']}"
                        new_record.append(entry)
                    if record['Organism']:
                        organism = f"organism={record['Organism']}"
                        new_record.append(organism)
                    if record['Transcript_id']:
                        transcript = f"transcript_id={record['Transcript_id']}"
                        new_record.append(transcript)
                    if record['Gene_id']:
                        gene_id = f"gene_id={record['Gene_id']}"
                        new_record.append(gene_id)

                    new_record = ";".join(new_record)
                    new_row = row.copy()
                    new_row[8] = new_record
                    w.write('\t'.join(new_row)+"\n")

    return bakta_gff_ext
