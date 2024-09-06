# %%
import sys
import os
import glob
import shutil
import subprocess
import pandas as pd
import re
import csv

# %%
input_path = sys.argv[1]
# input_path = "zvl_glu_ho_1"

# %%
os.getcwd()

# %%
# os.chdir(os.path.split(os.getcwd())[0])

# %%
# os.chdir("annotation_project")

# %%
for file in os.listdir(input_path):
    if file.endswith(".gff3"):
      bakta_gff = input_path + "/" + file
      bakta_tsv = input_path + "/" + file.rpartition('.')[0] + ".tsv"
      bakta_gff = os.path.abspath(bakta_gff)
      bakta_tsv = os.path.abspath(bakta_tsv)


# %%
if os.path.exists(bakta_gff):
    print(f'The file {bakta_gff} exists')
else:
    print(f'The file {bakta_gff} does not exist')
    exit()

# %%

bakta_gff_ext = bakta_gff.rpartition('.')[0] + "_extended.gff3"
bakta_tsv_ext = bakta_tsv.rpartition('.')[0] + "_extended.tsv"


# %%
shutil.copy(bakta_gff, bakta_gff_ext)
shutil.copy(bakta_tsv, bakta_tsv_ext)

# %%
pth = os.path.split(os.path.abspath(bakta_gff_ext))[0]

# %%
for file in os.listdir(pth):
    if file.endswith("ref2ref"):
        uni = os.path.abspath(pth) + "/" + file
for file in os.listdir(pth):
    if file.endswith("_uniref100_columns.tsv"):
        info_uniref100_table_tsv = os.path.abspath(pth) + "/" + file
        info_uniref100_table_ids = uni + "/uniprotinfo_uniref_representative_ids.tsv"
info_uniref100_tsv = uni + "/uniprotkb/uniprotinfo.tsv"
info_unknown_tsv = uni + "/uniprotkb/annotation/uniprotinfo.tsv"
info_unknown_table_tsv = uni + "/uniprotkb/annotation/UPIMAPI_results.tsv"

# %%
### adding unwkown proteins

# %%
unknown_df = pd.read_csv(info_unknown_table_tsv, sep = "\t", usecols=range(2), header = 0)

# %%
unknown_info = pd.read_csv(info_unknown_tsv, sep = "\t", header =0, index_col = 0)

# %%
pd.set_option('display.max_columns', None)

# %%
pd.set_option('display.max_rows', 30)

# %%
joined_unknown = unknown_df.join(unknown_info, on = "sseqid")

# %%
ext_tsv_df = pd.read_csv(bakta_tsv_ext, sep = "\t", header = None, comment = '#', \
                         names = ['Sequence Id','Type','Start','Stop','Strand','Locus Tag','Gene','Product','DbXrefs'])

# %%
joined_unknown['Gene Names'] = joined_unknown['Gene Names'].str.split().str[0]

# %%
merged_df = ext_tsv_df.merge(joined_unknown, left_on="Locus Tag", right_on='qseqid', suffixes=('', '_new'), how = "left")


# %%
joined_unknown.fillna('', inplace=True)

# %%
for i in range(len(merged_df)):
    # Get the current Locus Tag
    locus_tag = merged_df.at[i, 'Locus Tag']
    
    # Find the corresponding row in df2
    matching_row = joined_unknown[joined_unknown['qseqid'] == locus_tag]
    
    if not matching_row.empty and matching_row['Gene Names'].values[0] != '':
        # Update Gene Name and Product in df1
        merged_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
        merged_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
        merged_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
        merged_df.at[i, 'Entry UniProtKB'] = matching_row['sseqid'].values[0]

# %%
ext_tsv_df["Gene"] = merged_df["Gene"] 

# %%
ext_tsv_df["Product"] = merged_df["Product"] 

# %%
ext_tsv_df["Organism"] = merged_df["Organism"] 

# %%
ext_tsv_df.to_csv(bakta_tsv_ext, sep = '\t', index=False)

# %%
### adding uniref100 proteins

# %%
uniref_df_1 = pd.read_csv(info_uniref100_table_tsv, sep = "\t", header = None)

# %%
uniref_df_2 = pd.read_csv(info_uniref100_table_ids, sep = "\t", header = None)

# %%
uniref_df_2 = uniref_df_2.drop_duplicates(subset=0, keep='first')

# %%
joined_uni = uniref_df_1.merge(uniref_df_2, left_on =1, right_on = 0, how='left', suffixes=('', '_new'), validate = 'many_to_one')

# %%
joined_uni = joined_uni[["0", "1_new"]]

# %%
joined_uni.columns = ['id', 'uniprot']

# %%
uniref_df_3 = pd.read_csv(info_uniref100_tsv, sep = "\t", header = 0)

# %%
joined_uni = joined_uni.merge(uniref_df_3, left_on = "uniprot", right_on = "Entry", how='left', suffixes=('', '_new'))

# %%
joined_uni['Gene Names'] = joined_uni['Gene Names'].str.split().str[0]

# %%
ext_tsv_df = pd.read_csv(bakta_tsv_ext, sep = "\t", header = 0)

# %%
merged_df = ext_tsv_df.merge(joined_uni, left_on="Locus Tag", right_on='id', suffixes=('', '_new'), how = "left")


# %%
joined_uni.fillna('', inplace=True)

# %%
for i in range(len(merged_df)):
    # Get the current Locus Tag
    locus_tag = merged_df.at[i, 'Locus Tag']
    
    # Get UniProtKB from UserProtein
    uniprotkb = merged_df.at[i, 'DbXrefs']
    entry = re.compile(r"UserProtein:[^|]*\|([^,\n]*)")
    if isinstance(uniprotkb, str):
        match = re.search(entry, uniprotkb)
        if match:
            uniprotkb = match.group(1)
            merged_df.at[i, 'Entry UniProtKB'] = uniprotkb
            merged_df.at[i, 'Organism'] = "Escherichia coli"

    # Find the corresponding row in df2
    matching_row = joined_uni[joined_uni['id'] == locus_tag]

    if not matching_row.empty and matching_row['Gene Names'].values[0] != '':
        # Update Gene Name and Product in df1
        merged_df.at[i, 'Gene'] = matching_row['Gene Names'].values[0]
        merged_df.at[i, 'Product'] = matching_row['Protein names'].values[0]
        merged_df.at[i, 'Organism'] = matching_row['Organism'].values[0]
        merged_df.at[i, 'Entry UniProtKB'] = matching_row['Entry'].values[0]
    
    

# %%
ext_tsv_df["Gene"] = merged_df["Gene"] 

# %%
ext_tsv_df["Product"] = merged_df["Product"] 

# %%
ext_tsv_df["Organism"] = merged_df["Organism"] 

# %%
ext_tsv_df["Entry UniProtKB"] = merged_df["Entry UniProtKB"] 

# %%
ext_tsv_df.to_csv(bakta_tsv_ext, sep = '\t', index=False)

# %%
### changing GFF3

# %%
ext_tsv_df.fillna('', inplace=True)

# %%
with open(bakta_gff_ext, 'w') as w:
    with open(bakta_gff, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                w.write('\t'.join(row)+"\n")
                continue
            if len(row)<2:
                w.write('\t'.join(row)+"\n")
                continue              
            pairs = row[8].split(';')
            parsed_dict = {pair.split('=', 1)[0]: pair.split('=', 1)[1] for pair in pairs}
            id = ext_tsv_df[ ext_tsv_df['Locus Tag'] == parsed_dict.get('locus_tag') ]
            if id.empty:
                w.write('\t'.join(row)+"\n")
            else:
                record = id.iloc[0]
                new_record = []
                if record['Locus Tag']:
                    ID = f"ID={record['Locus Tag']}"
                    new_record.append(ID)
                if record['Product']:
                    Name = f"Name={record['Product'].replace(";", ",")}"
                    new_record.append(Name)
                if record['Locus Tag']:
                    locus_tag = f"locus_tag={record['Locus Tag']}"
                    new_record.append(locus_tag)
                if record['Product']:
                    product = f"product={record['Product'].replace(";", ",")}"
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
                
                new_record = ";".join(new_record)
                new_row = row.copy()
                new_row[8] = new_record
                w.write('\t'.join(new_row)+"\n")

                
