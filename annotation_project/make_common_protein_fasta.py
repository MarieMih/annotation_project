import os
import subprocess


def union_files(input_list, output_file):
    with open(output_file, 'w') as outfile:
        for fname in input_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def make_common_protein_fasta(tsvs):

    for i in ["_detected", "_unknown"]:
        fasta_files = list()
        for file in tsvs:
            fasta_files.append(file.replace("_extended.tsv", i + ".faa"))

        tmp = os.path.split(fasta_files[0])[0]
        tmp = os.path.split(tmp)[0]
        new_fasta_path = os.path.join(tmp, "union" + i + "_faa")
        if not os.path.exists(new_fasta_path):
            os.makedirs(new_fasta_path)
        new_fasta_file = os.path.join(new_fasta_path, "union" + i + ".faa")

        union_files(fasta_files, new_fasta_file)

        result = subprocess.run(['mmseqs', 'easy-cluster',
                                 new_fasta_file,
                                 os.path.join(new_fasta_path, "union"),
                                 os.path.join(new_fasta_path, "tmp"),
                                 "--cov-mode", "0",
                                 "-c", "0.8",
                                 "--min-seq-id", "0.9"],
                                check=True)

# integration test
# test_list = ["/storage/data1/marmi/annotation_project/rebuilding_data_test/bakta_annotation_GCF_000210475/GCF_000210475.tsv",
#              "/storage/data1/marmi/annotation_project/rebuilding_data_test/bakta_annotation_GCF_002853715/GCF_002853715.tsv",
#              "/storage/data1/marmi/annotation_project/rebuilding_data_test/bakta_annotation_GCF_003018455/GCF_003018455.tsv",
#              "/storage/data1/marmi/annotation_project/rebuilding_data_test/bakta_annotation_GCF_003697165/GCF_003697165.tsv"]

# make_common_protein_fasta(test_list)
