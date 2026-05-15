"""
Get all unique uniprotkb ids from tsvs and catch protein sequences.
"""
import os
import subprocess
import csv


def get_from_upimapi(id_set, directory, database, data_line, columns_to_return="Gene Names&Organism"):
    """
    Find information in db.
    """
    if database not in "UniProtKB,UniProtKB AC/ID,UniParc,UniRef50,UniRef90,UniRef100".split(","):
        print(f"{database} not valid UniProt database.")
        return

    tmpfile = "tmp_ids" + data_line + database.replace(" ","_") + ".csv"
    with open(tmpfile, "w") as outfile:
        outfile.write(",".join(id_set))
    subprocess.run(["upimapi",
                    "-i", "tmp_ids.csv",
                    "-o", directory,
                    "--columns", columns_to_return,
                    "--from-db", database],
                   check=True)
    os.remove(tmpfile)
