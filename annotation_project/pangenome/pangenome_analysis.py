import os
import itertools
import asyncio
import telegram_send
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from random import randint
import common_variables


def create_directory_with_soft_links(tsvs, target_or):
    target = os.path.abspath(target_or)  # отдебажить!!!
    if not os.path.exists(target):
        os.makedirs(target)
    for i in tsvs:
        new_link = os.path.split(i)[1]
        os.symlink(os.path.abspath(i), target + '/' + new_link)  # отдебажить!!!


def create_presence_absence_matrix(directory):
    gene_dict = {}

    for filename in os.listdir(directory):
        if filename.endswith('extended.tsv'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, sep='\t', header=None)
            gene_ids = df[6].unique()  # Gene symbol
            gene_dict[filename] = set(gene_ids)

    all_genes = list(set.union(*gene_dict.values()))
    presence_absence_matrix = pd.DataFrame(0, index=all_genes, columns=gene_dict.keys())

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
