import os
import csv
import pandas as pd
from find_unknown_proteins import finding_from_fasta


def create_acronym(phrase):
    if (phrase == "") or (phrase == "nan"):
        return "HP"
    trimmed = phrase.split('(')[0].strip()
    words = trimmed.split()
    acronym = "".join(word[0].upper() for word in words if word)
    return acronym


def correct_annotation_files(tsvs):
    """
    Get list of annotation tsv_extended-files and rewrite.
    """
    tmp = os.path.split(tsvs[0])[0]
    tmp = os.path.split(tmp)[0]
    fasta_file = os.path.join(tmp, "union_unknown_faa", "union_rep_seq.fasta")
    results = finding_from_fasta(fasta_file)

    cluster = fasta_file.replace("rep_seq.fasta", "cluster.tsv")
    cluster_df = pd.read_csv(cluster, sep="\t", header=None, names=['Head', 'Member'])
    cluster_df["File"] = cluster_df["Member"].str.rpartition("_")[0]

    for tsv in tsvs:
        gff = tsv.replace(".tsv", ".gff3")
        filename = os.path.split(gff)[-1].replace("_extended.gff3", "")
        cluster_part_df = cluster_df[cluster_df["File"].str.contains(filename, case=False)]
        df = pd.read_csv(tsv, sep="\t", header=0, comment='#')

        ### renaming gene_id and transcript_id

        df["Gene"] = df["Gene"].fillna('')
        df["Entry UniProtKB"] = df["Entry UniProtKB"].fillna('')
        df["Transcript_id"] = df["Type"] + "|" + df["Gene"].fillna('') + "|" + df["Entry UniProtKB"].fillna('')
        for i in range(len(df)):
            if ((df.at[i, 'Gene'] == "")
               and (df.at[i, 'Entry UniProtKB'] == "")
               and ((df.at[i, 'Type'] == "cds") or (df.at[i, 'Type'] == "sorf"))):
                df.at[i, 'Transcript_id'] = str(df.at[i, "Locus Tag"])
        df["Gene_id"] = df["Transcript_id"]

        ### rewrite tsv with finding information from uniprot

        for head in set(cluster_part_df["Head"]):

            if head in set(results["qseqid"]):
                r = results.loc[results["qseqid"] == head, :]
                for locustag in cluster_part_df.loc[cluster_part_df["Head"] == head, 'Member']:

                    if str(r["Protein names"].values[0]) != "nan":
                        df.loc[df['Locus Tag'] == locustag, ["Product"]] = r["Protein names"].values[0]

                    if str(r["Gene Names"].values[0]) != "nan":
                        df.loc[df['Locus Tag'] == locustag, ["Gene"]] = str(r["Gene Names"].values[0]).split(" ")[0]
                    else:
                        df.loc[df['Locus Tag'] == locustag, ["Gene"]] = "EGN_" + create_acronym(r["Protein names"].values[0]) + "_" + head

                    df.loc[df['Locus Tag'] == locustag, ["Organism"]] = r["Taxonomic lineage (SPECIES)"].values[0]
                    df.loc[df['Locus Tag'] == locustag, ["Entry UniProtKB"]] = r["sseqid"].values[0]
                    df.loc[df['Locus Tag'] == locustag, ['GO']] = r['Gene Ontology (GO)'].values[0]
                    df.loc[df['Locus Tag'] == locustag, ['KEGG']] = r['KEGG'].values[0]
                    df.loc[df['Locus Tag'] == locustag, ['UniPathway']] = r['UniPathway'].values[0]
                    df.loc[df['Locus Tag'] == locustag, ['Pathway']] = r['Pathway'].values[0]
                    df.loc[df['Locus Tag'] == locustag, ['Keywords']] = r['Keywords'].values[0]
                    df.loc[df['Locus Tag'] == locustag, ['Transcript_id']] = df.loc[df['Locus Tag'] == locustag, ["Type"]].values[0] + "|" + df.loc[df['Locus Tag'] == locustag, ["Gene"]].fillna('').values[0] + "|" + df.loc[df['Locus Tag'] == locustag, ["Entry UniProtKB"]].fillna('').values[0]
                    df.loc[df['Locus Tag'] == locustag, ['Gene_id']] = df.loc[df['Locus Tag'] == locustag, ['Transcript_id']].values[0]
            else:
                for locustag in cluster_part_df.loc[cluster_part_df["Head"] == head, 'Member']:
                    df.loc[df['Locus Tag'] == locustag, ['Transcript_id']] = head
                    df.loc[df['Locus Tag'] == locustag, ['Gene_id']] = head
                    df.loc[df['Locus Tag'] == locustag, ["Gene"]] = "UKW_HP_" + head

        ### writing into new tsv and gff

        new_tsv = tsv.replace(".tsv", "_polish.tsv")
        df.to_csv(new_tsv, sep='\t', index=False)
        df.fillna('', inplace=True)

        new_gff = gff.replace(".gff", "_polish.gff")
        with open(new_gff, 'w') as w:
            with open(gff.replace("_extended.gff3", ".gff3"), 'r') as f:

                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if row[0].startswith('#') or row[0].startswith('Sequence Id'):
                        w.write('\t'.join(row) + "\n")
                        continue
                    if len(row) < 2:
                        w.write('\t'.join(row) + "\n")
                        continue
                    pairs = row[8].split(';')
                    parsed_dict = {pair.split('=', 1)[0]: pair.split('=', 1)[1] for pair in pairs}
                    id = df[df['Locus Tag'] == parsed_dict.get('locus_tag')]
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
        os.remove(gff)
        os.rename(new_gff, gff)
        os.remove(tsv)
        os.rename(new_tsv, tsv)
