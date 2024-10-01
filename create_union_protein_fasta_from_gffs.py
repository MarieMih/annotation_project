"""
Get all unique uniprotkb ids from tsvs and catch protein sequences.
"""
import os
import subprocess
import csv


def get_from_upimapi(id_set, directory):
    """
    Find information in db.
    """
    with open("tmp_ids.csv", "w") as outfile:
        outfile.write(",".join(id_set))
    subprocess.run(["upimapi",
                    "-i", "tmp_ids.csv",
                    "-o", directory + "/upimapi_output",
                    "--fasta",
                    "--from-db", "UniProtKB AC/ID"],
                   check=True)
    os.remove("tmp_ids.csv")


def create_fasta_file(tsvs, directory):
    """
    Catch ids from all files and create fasta with upimapi
    """
    ids = set()
    for i in tsvs:
        with open(i, "r", newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)
            for row in reader:
                ids.add(row[10])
    ids.remove("")
    get_from_upimapi(ids, directory)
