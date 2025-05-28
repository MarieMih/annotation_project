import os
import csv
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

COLUMNS = ["key", "sample_name", 'empty_gene_name', 'non_empty_gene_name', 'gene_name_uni', 'user_protein', 'total_lines', 'uncharacterized_protein', 'hypothetical_protein', 'uniprotkb', 'uniprotkb_uni']
# COLUMNS = ["key", "sample_name", 'empty_gene_name', 'non_empty_gene_name', 'gene_name_uni', 'user_protein', 'total_lines', 'uncharacterized_protein', 'hypothetical_protein']


def process_tsv(input_file):
    counts = defaultdict(lambda: {
        'sample_name': os.path.splitext(os.path.split(input_file)[1])[0],
        'empty_sixth': 0,
        'non_empty_sixth': 0,
        'gene_uni': set(),
        'user_protein': 0,
        'total_lines': 0,
        'uncharacterized_protein': 0,
        'hypothetical_protein': 0,
        'uniprotkb': 0,
        'uniprotkb_uni': set()
    })

    df = pd.DataFrame(columns=COLUMNS)
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')

        for row in reader:
            if row[0].startswith('#'):
                continue

            key = row[1]
            counts[key]['total_lines'] += 1

            if row[6] == '':
                counts[key]['empty_sixth'] += 1
            else:
                counts[key]['non_empty_sixth'] += 1
                counts[key]['gene_uni'].add(row[6])

                if "UserProtein" in row[8]:
                    counts[key]['user_protein'] += 1

                if "Uncharacterized protein" in row[7]:
                    counts[key]['uncharacterized_protein'] += 1

            if "hypothetical protein" in row[7]:
                counts[key]['hypothetical_protein'] += 1

            if row[10] != '':
                counts[key]['uniprotkb'] += 1
                counts[key]['uniprotkb_uni'].add(row[10])

    print(f"{os.path.split(input_file)[1]}")

    for key, count in sorted(counts.items(), key=lambda item: item[0]):
        # print(f"Group {key}:")
        # print(f"  Total lines: {count['total_lines']}")
        # print(f"  Empty gene names: {count['empty_sixth']}")
        # print(f"  Non-empty gene names: {count['non_empty_sixth']}")
        # print(f"  Non-empty unique gene names: {len(count['gene_uni'])}")
        # print(f"  Non-empty gene names and 'UserProtein': {count['user_protein']}")
        # print(f"  Non-empty gene names and 'Uncharacterized protein': {count['uncharacterized_protein']}")
        # print(f"  Lines with 'hypothetical protein': {count['hypothetical_protein']}")
        # print(f"  UniProtKB: {count['uniprotkb']}")
        # print(f"  UniProtKB (unique): {len(count['uniprotkb_uni'])}")
        # print()
        df.loc[len(df.index)] = [key, count['sample_name'], count['empty_sixth'], count['non_empty_sixth'], len(count['gene_uni']), count['user_protein'], count['total_lines'], count['uncharacterized_protein'], count['hypothetical_protein'], count['uniprotkb'], len(count['uniprotkb_uni'])]
        # df.loc[len(df.index)] = [key, count['sample_name'], count['empty_sixth'], count['non_empty_sixth'], len(count['gene_uni']), count['user_protein'], count['total_lines'], count['uncharacterized_protein'], count['hypothetical_protein']]
    df = df[df.key != "Type"]

    return df

def make_stat_file(directory, filename=None):
    files = []

    for curfile in os.listdir(directory):
        if curfile.endswith('.tsv') and not curfile.startswith('stat_annotation'):
            file_path = os.path.join(directory, curfile)
            files.append(file_path)

    df = pd.DataFrame([], columns=COLUMNS)
    for i in files:
        df = pd.concat([df, process_tsv(i)])

    if filename is None:
        file_output = os.path.join(directory, "stat_annotation.tsv")
    else:
        file_output = filename
    df.to_csv(file_output, sep='\t', index=False)

    df = df[df.key == "cds"]
    df = df.drop("key", axis=1)

    df.set_index('sample_name', inplace=True)
    df = df.apply(pd.to_numeric, errors='coerce')

    plt.figure(figsize=(8, 6))
    sns.heatmap(df, annot=True, fmt="g", cmap="viridis", cbar_kws={'label': 'Gene Count'}, linewidths=0.5)

    plt.title('CDS Gene Count Heatmap')
    plt.xlabel('Gene Type')
    plt.ylabel('Sample Names')
    plt.xticks(rotation=20, ha='right')
    plt.tight_layout()

    plt.savefig(os.path.join(directory, "stat_annotation.png"), bbox_inches='tight')
