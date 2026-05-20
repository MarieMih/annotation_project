import os
import subprocess
import common_variables
from helpers import union_files
import pandas as pd
from Bio import SeqIO


def make_feature_clusters(tool, tsv, dir, outdir, fname="rna_clusters"):
    """Заглушка для некодирующих фичей"""

    new_fasta_file = os.path.join(dir, "union_noncds.fna")
    ids = make_noncoding_rna_fasta(tsv)
    fna = [i.replace(".tsv", ".rna.ffn") for i in tsv]
    union_files(fna, new_fasta_file)
    df = pd.DataFrame({
        'parent': ids,
        'child': ids
    })
    
    cluster_file = os.path.join(outdir, f"{fname}.tsv")
    df.to_csv(cluster_file, header=False, sep="\t", index=False)
    return cluster_file


def make_noncoding_rna_fasta(tsv, type="crispr crispr-repeat crispr-spacer ncRNA ncRNA-region oriC oriT oriV rRNA tRNA tmRNA", feature="rna"):
    names = common_variables.BAKTA_TSV_HEADER.split(",")

    out_rec = []
    common_ids = []
    for file in tsv:
        df = pd.read_csv(file, sep="\t", header=0, comment="#", names=names)
        df = df[df["Type"].isin(type.split())]
        ids = list(df["Locus Tag"].values)
        with open(file.replace(".tsv", ".ffn")) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in ids:
                    out_rec.append(record)
        common_ids.extend(ids)
        with open(file.replace(".tsv", f".{feature}.ffn"), "w") as output_handle:
            SeqIO.write(out_rec, output_handle, "fasta")
    return(common_ids)
