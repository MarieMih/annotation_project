import os
import itertools
import asyncio
import telegram_send
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint
import common_variables
from collections import defaultdict


def create_presence_absence_matrix(directory, cluster_file, mode, output_directory):
    names = "Sequence Id,Type,Start,Stop,Strand,Locus Tag,Gene,Product,DbXrefs".split(",")
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.tsv'):
            files.append(filename.replace(".tsv",""))

    clusters  = pd.read_csv(cluster_file, sep='\t', header=None, names=["parent", "child"])
    parent_locus = list(set(clusters["parent"]))
    clusters  = dict(zip(clusters["child"], clusters["parent"]))

    match mode:
        case "binary":
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
            df = pd.read_csv(file_path, sep='\t', comment="#", header=None, names=names)
            df = df[df["Type"].isin("cds sorf".split())]
            gene_ids = set(df["Locus Tag"].unique())

            if mode == "locus":
                parent_fasta[filename].update(set.intersection(set(parent_locus), gene_ids))

            for i in gene_ids:
                match mode:
                    case "binary":
                        value = 1
                    case "locus":
                        value = i
                presence_absence_matrix.loc[clusters[i], filename.replace(".tsv", "")] = value 

    match mode:
        case "binary":
            presence_absence_matrix['count'] = presence_absence_matrix.sum(axis=1)
        case "locus":
            presence_absence_matrix['count'] = (presence_absence_matrix != "").sum(axis=1)
    
    presence_absence_matrix = presence_absence_matrix.sort_values(by='count', ascending=False)
    presence_absence_matrix = presence_absence_matrix.drop(columns=['count'])

    if mode == "binary":
        presence_absence_matrix = presence_absence_matrix.reset_index(names='Locus Tag')
        presence_absence_matrix.to_csv(os.path.join(output_directory, 'presence_absence_matrix_binary.tsv'), index=False, sep="\t")
        return os.path.join(output_directory, 'presence_absence_matrix_binary.tsv')

    if mode == "locus":
        presence_absence_matrix["Gene"] = ""
        presence_absence_matrix["Product"] = ""
    
        for i in parent_fasta.keys():
            file_path = os.path.join(directory, i)
            df = pd.read_csv(file_path, sep='\t', comment="#", header=None, names=names)  
            to_update = df[df["Locus Tag"].isin(parent_fasta[i])]
            to_update = to_update.set_index("Locus Tag")
            to_update = to_update["Gene Product".split()]
            presence_absence_matrix.loc[to_update.index, ['Gene', 'Product']] = to_update[['Gene', 'Product']]
            presence_absence_matrix["Genome"] = i.replace(".tsv", "")

        presence_absence_matrix = presence_absence_matrix.reset_index(names='Locus Tag')
        presence_absence_matrix["PID"] = range(1, len(presence_absence_matrix) + 1)

        cols_to_move = ['Locus Tag', "Genome", 'PID', "Gene", "Product"]
        new_order = cols_to_move + [c for c in presence_absence_matrix.columns if c not in cols_to_move]
        presence_absence_matrix = presence_absence_matrix[new_order]

        presence_absence_matrix.to_csv(os.path.join(output_directory, 'presence_absence_matrix.tsv'), index=False, sep="\t")
        return os.path.join(output_directory, 'presence_absence_matrix.tsv')      



def create_presence_absence_matrix_by_symbol(directory):
    gene_dict = {}

    for filename in os.listdir(directory):
        if filename.endswith('extended.tsv'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, sep='\t', header=None)
            gene_ids = df["Gene"].unique()  # Gene symbol
            gene_dict[filename] = set(gene_ids)

    parent_locus = list(set.union(*gene_dict.values()))
    presence_absence_matrix = pd.DataFrame(0, index=parent_locus, columns=gene_dict.keys())

    for filename, genes in gene_dict.items():
        for i in genes:
            presence_absence_matrix.loc[i, filename] = 1

    presence_absence_matrix['count'] = presence_absence_matrix.sum(axis=1)
    presence_absence_matrix = presence_absence_matrix.sort_values(by='count', ascending=False)
    presence_absence_matrix = presence_absence_matrix.drop(columns=['count'])
    presence_absence_matrix.to_csv(directory + "/" + 'presence_absence_matrix.csv')

    return directory + "/" + 'presence_absence_matrix.csv'


def calculate_core_genome_combinations(df):
    core_genome_sizes = {i: [] for i in range(1, len(df.columns) + 1)}

    for i in range(1, len(df.columns) + 1):
        count = 0
        for combo in itertools.combinations(df.columns, i):
            subset = df[list(combo)]
            core_genes = subset.sum(axis=1) == i
            core_genome_sizes[i].append(core_genes.sum())
            if count == 1000:
                break
            count += 1

    return core_genome_sizes


def calculate_pangenome_combinations(df):
    pangenome_sizes = {i: [] for i in range(1, len(df.columns) + 1)}

    for i in range(1, len(df.columns) + 1):
        count = 0
        for combo in itertools.combinations(df.columns, i):
            subset = df[list(combo)]
            pangenome_genes = subset.sum(axis=1) > 0
            pangenome_sizes[i].append(pangenome_genes.sum())
            if count == 1000:
                break
            count += 1

    return pangenome_sizes


async def send_smth(cor_image, pan_image):
    with open(cor_image, "rb") as f:
        await telegram_send.send(images=[f])
    with open(pan_image, "rb") as f:
        await telegram_send.send(images=[f])


def pangenome_analysis(directory_or):
    directory = os.path.abspath(directory_or)  # отдебажить!!!
    file_path = create_presence_absence_matrix(directory)
    data = pd.read_csv(file_path, sep=',', index_col=0)
    data.columns = [i.replace("_extended", "") for i in data.columns]

    core_genome_sizes = calculate_core_genome_combinations(data)
    boxplot_data = [core_genome_sizes[i] for i in range(1, len(data.columns) + 1)]

    plt.figure(figsize=(10, 6))
    plt.boxplot(boxplot_data, tick_labels=[str(i) for i in range(1, len(data.columns) + 1)])
    plt.xlabel('Number of Samples')
    plt.ylabel('Core Genome Size')
    plt.title('Core Genome Size Distribution by Number of Samples')
    plt.grid(True)
    cor_image = directory + "/" + "cor.png"
    plt.savefig(cor_image)

    pangenome_sizes = calculate_pangenome_combinations(data)
    boxplot_data = [pangenome_sizes[i] for i in range(1, len(data.columns) + 1)]

    plt.figure(figsize=(10, 6))
    plt.boxplot(boxplot_data, tick_labels=[str(i) for i in range(1, len(data.columns) + 1)])
    plt.xlabel('Number of Samples')
    plt.ylabel('Pangenome Size')
    plt.title('Pangenome Size Distribution by Number of Samples')
    plt.grid(True)
    pan_image = directory + "/" + "pan.png"
    plt.savefig(pan_image)

    if common_variables.SEND_NOTIFICATION:
        asyncio.run(send_smth(cor_image, pan_image))



def pangenome_tsv(directory, cluster_file, matrix, output_directory):
    names = "Sequence Id,Type,Start,Stop,Strand,Locus Tag,Gene,Product,DbXrefs".split(",")
    pannames = "Sequence Id,Type,Start,Stop,Strand,Genome,Gene,Gene Name,Gene synonymes,Product,DbXrefs,Organism,Inference,KEGG,GO,Gene id,Transcript id".split(",")

    clusters  = pd.read_csv(cluster_file, sep='\t', header=None, names=["parent", "child"])
    clusters  = dict(zip(clusters["child"], clusters["parent"]))

    pangenome_table = pd.read_csv(matrix, sep="\t", index_col=0, header=0)
    pangenome_table = pangenome_table[["PID", "Genome"]]
    pangenome_table[pannames] = ""

    parent_fasta = pangenome_table.groupby('Genome').apply(lambda x: set(x.index)).to_dict()
    
    for i in parent_fasta.keys():
        file_path = os.path.join(directory, i + ".tsv")
        df = pd.read_csv(file_path, sep='\t', comment="#", header=None, names=names)  
        to_update = df[df["Locus Tag"].isin(parent_fasta[i])]
        to_update = to_update.set_index("Locus Tag")
        pangenome_table.loc[to_update.index, to_update.columns] = to_update

    pangenome_table = pangenome_table.reset_index(names='PID Locus Tag')

    # cols_to_move = ['Locus Tag', 'PID', "Gene", "Product"]
    # new_order = cols_to_move + [c for c in pangenome_table.columns if c not in cols_to_move]
    # pangenome_table = pangenome_table[new_order]

    pangenome_table.to_csv(os.path.join(output_directory, 'pangenome_table.tsv'), index=False, sep="\t")
    return os.path.join(output_directory, 'pangenome_table.tsv')   
