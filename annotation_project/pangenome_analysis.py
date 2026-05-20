import os
import itertools
import asyncio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from math import comb
from random import randint, sample, shuffle
import common_variables
from collections import defaultdict
from helpers import create_acronym, send_smth
from Bio import SeqIO


def create_presence_absence_matrix(directory, cluster_file, mode, output_directory, fname = 'presence_absence_matrix'):

    names = common_variables.BAKTA_TSV_HEADER.split(",")

    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.tsv'):
            files.append(filename.replace(".tsv",""))

    clusters     = pd.read_csv(cluster_file, sep='\t', header=None, names=["parent", "child"])
    parent_locus = list(set(clusters["parent"]))
    clusters     = dict(zip(clusters["child"], clusters["parent"]))

    match mode:
        case "binary":
            initial_value = 0
        case "numeric":
            initial_value = 0
        case "locus":
            initial_value = ""
            parent_fasta = defaultdict(set)
        case _:
            print(f"Not correct mode.")
            return

    presence_absence_matrix = pd.DataFrame(initial_value, index=parent_locus, columns=files)

    for filename in os.listdir(directory):
        if filename.endswith('.tsv'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, sep='\t', comment="#", header=0, names=names)
            gene_ids = set(df["Locus Tag"].unique())

            if mode == "locus":
                parent_fasta[filename].update(set.intersection(set(parent_locus), gene_ids))

            for i in gene_ids:
                if clusters.get(i) is not None:
                    match mode:
                        case "binary":
                            value = 1
                            presence_absence_matrix.loc[clusters[i], filename.replace(".tsv", "")] = value
                        case "numeric":
                            value = 1
                            presence_absence_matrix.loc[clusters[i], filename.replace(".tsv", "")] += 1
                        case "locus":
                            value = i
                            current_value = presence_absence_matrix.loc[clusters[i], filename.replace(".tsv", "")]
                            if current_value == "":
                                presence_absence_matrix.loc[clusters[i], filename.replace(".tsv", "")] = value
                            else:
                                presence_absence_matrix.loc[clusters[i], filename.replace(".tsv", "")] = current_value + ";" + value

    match mode:
        case "binary":
            presence_absence_matrix['count'] = presence_absence_matrix.sum(axis=1)
        case "numeric":
            presence_absence_matrix['count'] = presence_absence_matrix.sum(axis=1)
        case "locus":
            presence_absence_matrix['count'] = (presence_absence_matrix != "").sum(axis=1)
    
    presence_absence_matrix = presence_absence_matrix.sort_values(by='count', ascending=False)
    presence_absence_matrix = presence_absence_matrix.drop(columns=['count'])

    if mode == "binary":
        presence_absence_matrix = presence_absence_matrix.reset_index(names='Locus Tag')
        presence_absence_matrix.to_csv(os.path.join(output_directory, fname + '_binary.tsv'), index=False, sep="\t")

    if mode == "numeric":
        presence_absence_matrix = presence_absence_matrix.reset_index(names='Locus Tag')
        presence_absence_matrix.to_csv(os.path.join(output_directory, fname + '_numeric.tsv'), index=False, sep="\t")

    if mode == "locus":
        presence_absence_matrix["Gene"] = ""
        presence_absence_matrix["Product"] = ""
    
        for i in parent_fasta.keys():
            file_path = os.path.join(directory, i)
            df = pd.read_csv(file_path, sep='\t', comment="#", header=0, names=names)  
            to_update = df[df["Locus Tag"].isin(parent_fasta[i])]
            to_update = to_update.set_index("Locus Tag")
            to_update = to_update["Gene Product".split()]
            presence_absence_matrix.loc[to_update.index, ['Gene', 'Product']] = to_update[['Gene', 'Product']]
            presence_absence_matrix.loc[to_update.index, ['Genome']] = i.replace(".tsv", "")

        presence_absence_matrix = presence_absence_matrix.reset_index(names='Locus Tag')
        presence_absence_matrix["PID"] = range(1, len(presence_absence_matrix) + 1)

        cols_to_move = ['Locus Tag', "Genome", 'PID', "Gene", "Product"]
        new_order = cols_to_move + [c for c in presence_absence_matrix.columns if c not in cols_to_move]
        presence_absence_matrix = presence_absence_matrix[new_order]

        presence_absence_matrix.to_csv(os.path.join(output_directory, fname + '.tsv'), index=False, sep="\t") 

    match mode:
        case "binary":
            return os.path.join(output_directory, fname + '_binary.tsv')
        case "numeric":
            return os.path.join(output_directory, fname + '_numeric.tsv')
        case "locus":
            return os.path.join(output_directory, fname + '.tsv')  


def pangenome_tsv(directory, cluster_file, matrix, output_directory, fname = 'pangenome_table'):
    names      = common_variables.BAKTA_TSV_HEADER.split(",")
    infr_names = common_variables.BAKTA_INFERENCE_TSV_HEADER.split(",")
    pannames   = common_variables.PANGENOME_TSV_HEADER.split(",")

    clusters = pd.read_csv(cluster_file, sep='\t', header=None, names=["parent", "child"])
    clusters = dict(zip(clusters["child"], clusters["parent"]))

    pangenome_table = pd.read_csv(matrix, sep="\t", index_col=0, header=0)
    pangenome_table = pangenome_table[["PID", "Genome"]]
    pangenome_table[pannames] = ""
    pangenome_table[["Start", "Stop"]] = -1
    pangenome_table = pangenome_table.astype({"Start": int, "Stop": int})

    parent_fasta = pangenome_table.groupby('Genome').apply(lambda x: set(x.index)).to_dict()
    
    for i in parent_fasta.keys():
        file_path = os.path.join(directory, i + ".tsv")
        infr_path = os.path.join(directory.replace("tsvs", "inferences"), i + ".inference.tsv")

        df = pd.read_csv(file_path, sep='\t', comment="#", header=0, names=names)  
        to_update = df[df["Locus Tag"].isin(parent_fasta[i])]
        to_update = to_update.set_index("Locus Tag")
        pangenome_table.update(to_update)

        df = pd.read_csv(infr_path, sep='\t', comment="#", header=0, names=infr_names)  
        to_update = df[df["Locus Tag"].isin(parent_fasta[i])]
        to_update = to_update.set_index("Locus Tag")
        to_update = to_update.rename(columns={"Accession": "Inference"})
        to_update = to_update["Inference"]
        pangenome_table.update(to_update)

    # check_api = pangenome_table["Inference"].str.split(":")

    pangenome_table["Gene Name"] = np.where(
        pangenome_table["Gene"] == "",
        "EGN_" + pangenome_table["Product"].apply(create_acronym) + "_P" + pangenome_table["PID"].astype("str"),
        pangenome_table["Gene"]
    )

    pangenome_table["gene_id"]       = pangenome_table["Type"] + "|" + pangenome_table["Gene Name"] + "|" + pangenome_table["Inference"].astype("str")
    pangenome_table["transcript_id"] = pangenome_table["gene_id"]

    pangenome_table = pangenome_table.reset_index(names='PID Locus Tag')

    pangenome_table.to_csv(os.path.join(output_directory, fname + '.tsv'), index=False, sep="\t")
    return os.path.join(output_directory, fname + '.tsv')   


def pangenome_fasta(directory, pangenome_table, output_directory, ftype: str, fname = "pangenome"):
    if ftype not in "faa ffn".split():
        print(f"{ftype} not fasta.")
        return

    pangenome_table = pd.read_csv(pangenome_table, sep="\t", index_col=None, header=0)
    parent_fasta = pangenome_table.groupby('Genome')["PID Locus Tag"].apply(set).to_dict()

    sorted_tags = list(pangenome_table["PID Locus Tag"])

    out_rec = defaultdict()
    for i in parent_fasta.keys():
        tags    = parent_fasta[i]
        with open(os.path.join(directory, i + "." + ftype)) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in tags:
                    out_rec[record.id] = record

    sorted_data = dict(sorted(out_rec.items(), key=lambda x: sorted_tags.index(x[0])))

    with open(os.path.join(output_directory, fname + "." + ftype), "w") as output_handle:
        SeqIO.write(list(sorted_data.values()), output_handle, "fasta")


def calculate_core_genome_combinations(df):
    core_genome_sizes = {i: [] for i in range(1, len(df.columns) + 1)}

    for i in range(1, len(df.columns) + 1):
        if comb(len(df.columns), i) < 1000:
            for combo in itertools.combinations(df.columns, i):
                subset = df[list(combo)]
                core_genes = subset.sum(axis=1) == i
                core_genome_sizes[i].append(core_genes.sum())
        else:
            count = 0
            for _ in range(1000):
                subset = df[list(sample(sorted(df.columns), i))]
                core_genes = subset.sum(axis=1) == i
                core_genome_sizes[i].append(core_genes.sum())
                if count == 1000:
                    break
                count += 1

    return core_genome_sizes


def calculate_pangenome_combinations(df):
    pangenome_sizes = {i: [] for i in range(1, len(df.columns) + 1)}

    for i in range(1, len(df.columns) + 1):
        if comb(len(df.columns), i) < 1000:
            for combo in itertools.combinations(df.columns, i):
                subset = df[list(combo)]
                pangenome_genes = subset.sum(axis=1) > 0
                pangenome_sizes[i].append(pangenome_genes.sum())
        else:
            count = 0
            for _ in range(min(comb(len(df.columns), i), 1000)):
                subset = df[list(sample(sorted(df.columns), i))]
                pangenome_genes = subset.sum(axis=1) > 0
                pangenome_sizes[i].append(pangenome_genes.sum())
                if count == 1000:
                    break
                count += 1

    return pangenome_sizes





def pangenome_curves(file_path, corefname="core_genome_size_distribution", panfname="pangenome_size_distribution"):
    directory = os.path.dirname(file_path)
    data = pd.read_csv(file_path, sep='\t', index_col=0, header=0)
    set_size = len(data.columns)

    core_genome_sizes = calculate_core_genome_combinations(data)

    boxplot_data = [core_genome_sizes[i] for i in range(1, set_size + 1)]
    plt.figure(figsize=(10, 6))
    # sns.boxplot(boxplot_data, color="0.8", width=0.3)
    # sns.stripplot(boxplot_data)
    sns.violinplot(boxplot_data)
    plt.xlabel('Number of Samples')
    plt.ylabel('Core Genome Size')
    plt.title('Core Genome Size Distribution by Number of Samples')
    plt.grid(True)
    plt.xticks(
        ticks=range(set_size),
        labels=range(1, set_size + 1)
    )
    cor_image = os.path.join(directory, corefname + '.png')
    plt.savefig(cor_image)

    pangenome_sizes = calculate_pangenome_combinations(data)

    boxplot_data = [pangenome_sizes[i] for i in range(1, set_size + 1)]
    plt.figure(figsize=(10, 6))
    # sns.boxplot(boxplot_data, color="0.8", width=0.3)
    # sns.stripplot(boxplot_data)
    sns.violinplot(boxplot_data)
    plt.xlabel('Number of Samples')
    plt.ylabel('Pangenome Size')
    plt.title('Pangenome Size Distribution by Number of Samples')
    plt.grid(True)
    plt.xticks(
        ticks=range(set_size),
        labels=range(1, set_size + 1)
    )
    pan_image = os.path.join(directory, panfname + '.png')
    plt.savefig(pan_image)

    if common_variables.SEND_NOTIFICATION:
        asyncio.run(send_smth(cor_image, pan_image))

    return [core_genome_sizes[set_size][0], pangenome_sizes[set_size][0]]
