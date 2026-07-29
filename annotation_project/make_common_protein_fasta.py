import os
import subprocess
import shutil
import csv
from helpers import union_files


def make_cds_clusters(tool, faa, dir, outdir, fname="clusters", path_to_gff=None):
    match tool:
        case "MMSEQS2":
            cluster_file = make_with_MMSEQS2(faa, dir)
        case "PGAP2":
            cluster_file = make_with_PGAP2(path_to_gff, dir)
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

def make_with_PGAP2(path_to_gff, dir):
    new_path = os.path.join(dir, "output_pgap2")
    if not os.path.exists(new_path):
        os.makedirs(new_path)

    result = subprocess.run(['pgap2', 'main',
                                "-i", path_to_gff,
                                "-o", new_path],
                            check=True)
    
    parse_clusters(os.path.join(new_path, "pgap2.partition.gene_content.csv"),
                   os.path.join(new_path, "pgap2.partition.gene_content.parsed.tsv"))

    return os.path.join(new_path, "pgap2.partition.gene_content.parsed.tsv")

def parse_clusters(input_file, output_file):
    with open(input_file, newline="") as f_in, open(output_file, "w", newline="") as f_out:
        reader = csv.reader(f_in)
        writer = csv.writer(f_out, delimiter="\t")

        next(reader)

        for row in reader:
            loci = []

            for value in row[1:]:
                if not value:
                    continue

                for locus in value.split(";"):
                    locus = locus.strip()
                    if locus:
                        loci.append(locus)

            if not loci:
                continue

            representative = loci[0]

            for member in loci:
                writer.writerow([representative, member])
