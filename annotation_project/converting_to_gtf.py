import os
import subprocess
import pandas as pd
import re


def convert_AGAT(file):
    """
    Use AGAT.
    """
    new_gtf = file.rpartition(".")[0] + ".gtf"
    subprocess.run(["agat_convert_sp_gff2gtf.pl",
                    "--gff", os.path.abspath(file),  # отдебажить!!
                    "-o", os.path.abspath(new_gtf)],
                   check=True,
                   cwd=os.path.split(file)[0])


def convert_gff_to_gtf(file):

    convert_AGAT(file)
    new_gtf = file.rpartition(".")[0] + ".gtf"

    commented_lines = []
    with open(new_gtf, "r") as file:
        for line in file:
            if line.startswith("#"):
                commented_lines.append(line)
            else:
                break

    col_names = ["seqid", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    gtf_data = pd.read_csv(new_gtf, sep="\t", comment="#", header=None, names=col_names)

    def clean_multiple_quotes(attributes):
        fields = attributes.split(";")
        cleaned_fields = []

        for field in fields:
            if field.count('"') > 2:
                first_quote_match = re.search(r'"[^"]+"', field)
                if first_quote_match:
                    cleaned_field = field[:first_quote_match.end()] + re.sub(r'"[^"]+"', '', field[first_quote_match.end():])
                    cleaned_fields.append(cleaned_field.strip())
            else:
                cleaned_fields.append(field.strip())

        return "; ".join(cleaned_fields).strip()

    gtf_data['attributes'] = gtf_data['attributes'].apply(clean_multiple_quotes)
    gtf_data.loc[gtf_data['feature'] == 'CRISPR', 'strand'] = '+'

    output_gtf = new_gtf.replace(".gtf", "_tmp.gtf")
    with open(output_gtf, "w") as file:
        file.writelines(commented_lines)
        gtf_data.to_csv(file, sep="\t", index=False, header=False, quoting=3)

    os.remove(new_gtf)
    os.rename(output_gtf, new_gtf)
