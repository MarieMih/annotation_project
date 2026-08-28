import os
import csv
import pandas as pd
import numpy as np
from helpers import create_acronym, check_file_exists

def precorrect_tsv_file(tsv_path):

    check_file_exists(tsv_path)

    names = "Sequence Id,Type,Start,Stop,Strand,Locus Tag,Gene,Product,DbXrefs".split(",")

    # if check_file_exists(tsv_path.replace(".tsv", ".bakta.tsv")) != 0:
    #     origin = pd.read_csv(tsv_path, sep="\t", comment="#", names=names, header=None)
    #     print(f"{tsv_path} 1")
    # else:
    #     origin = pd.read_csv(tsv_path, sep="\t", comment="#", names=names, header=0)
    #     print(f"{tsv_path} 2")

    with open(tsv_path, "r") as f:
        first_line = f.readline()

    if first_line.startswith("#"):
        origin = pd.read_csv(
            tsv_path,
            sep="\t",
            comment="#",
            names=names,
            header=None
        )
    else:
        origin = pd.read_csv(
            tsv_path,
            sep="\t",
            comment="#",
            names=names,
            header=0
        )

    def count_calls(func):
        def wrapper(*args, **kwargs):
            wrapper.count += 1
            return func(*args, **kwargs)
        wrapper.count = 0
        return wrapper
    
    @count_calls
    def new_feature():
        return

    def feature_name(genome, prefix = "NONCDSF"):
        new_feature()
        return f"{prefix}_{genome[-16:]}_{new_feature.count}"
    
    mask = ((origin["Locus Tag"] == "") | (origin["Locus Tag"].isna()))
    origin.loc[mask, "Locus Tag"] = (
        origin.loc[mask, "Locus Tag"]
        .map(lambda _: feature_name(os.path.basename(tsv_path).replace(".tsv", "")))
    )

    origin.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


def correct_tsv_file(tsv_path, pangenome_table, cluster_file):

    check_file_exists(tsv_path)

    if "_pangenome.tsv" not in tsv_path:
        output_tsv = tsv_path.replace(".tsv", "_pangenome.tsv")
    else:
        output_tsv = tsv_path

    origin = pd.read_csv(tsv_path, sep="\t", comment="#", header=0, dtype={"PID": str})
    if "Uniq Gene" not in origin.columns:
        origin = origin.rename(columns={"Gene": "Uniq Gene", "Product": "Uniq Product", "DbXrefs": "Uniq DbXrefs"})

    pangenome_table_df = pd.read_csv(pangenome_table, sep="\t", index_col=None, header=0, dtype={"PID": str})

    clusters  = pd.read_csv(cluster_file, sep='\t', header=None, names=["parent", "child"])
    clusters  = dict(zip(clusters["child"], clusters["parent"]))

    origin["Parent LT"] = ""
    origin["Parent LT"] = origin["Locus Tag"].map(clusters).fillna(origin["Parent LT"])

    pangenome_table_df = pangenome_table_df.drop("Sequence Id,Type,Start,Stop,Strand".split(","), axis=1)

    if "PID Locus Tag" not in origin.columns:
        origin = pd.merge(origin, pangenome_table_df, left_on='Parent LT', right_on='PID Locus Tag', how='left')
        origin = origin.drop(columns =["Parent LT"])
    else:
        tmp_df = pd.merge(origin, pangenome_table_df, left_on='Parent LT', right_on='PID Locus Tag', how='left', suffixes=["_old", ""], sort=False)
        tmp_df.index = origin.index
        tmp_df = tmp_df[tmp_df["PID Locus Tag"] != ""]
        cols_to_update = [
            col for col in pangenome_table_df.columns
            if col not in ["PID Locus Tag", "PID"] ########################## when PID of different features will be unique, "PID" could be erased
        ]
        origin.update(tmp_df[cols_to_update])
        origin = origin.drop(columns =["Parent LT"])


    cols_to_move = ['DbXrefs', "Uniq DbXrefs", "PID Locus Tag", "PID", "Genome"]
    new_order = [c for c in origin.columns if c not in cols_to_move] + cols_to_move
    origin = origin[new_order]

    origin.to_csv(output_tsv, sep="\t", index=False)
    return output_tsv
    

