import os
import subprocess
import shutil
from helpers import union_files


def make_cds_clusters(tool, faa, dir, outdir, fname="clusters"):
    match tool:
        case "MMSEQS2":
            cluster_file = make_with_MMSEQS2(faa, dir)
        case _:
            print(f"{tool} is not correct value and not supported now.")
            return None
    final_cluster_file = os.path.join(outdir, f"{fname}.tsv")
    shutil.copy(cluster_file, final_cluster_file)
    return final_cluster_file


def make_with_MMSEQS2(fasta_files, dir):
    new_fasta_path = os.path.join(dir, "mmseqs_results_faa")
    if not os.path.exists(new_fasta_path):
        os.makedirs(new_fasta_path)
    new_fasta_file = os.path.join(new_fasta_path, "union.faa")

    union_files(fasta_files, new_fasta_file)

    result = subprocess.run(['mmseqs', 'easy-cluster',
                                new_fasta_file,
                                os.path.join(new_fasta_path, "union"),
                                os.path.join(new_fasta_path, "tmp"),
                                "--cov-mode", "0",
                                "-c", "0.8",
                                "--min-seq-id", "0.9"],
                            check=True)
    
    return os.path.join(new_fasta_path, "union_cluster.tsv")
