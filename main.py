import os
import sys
import asyncio
from preparation import bakta_annotation
from preparation import send_smth
from annotation import annotation
from pangenome.pangenome_analysis import create_directory_with_soft_links
from pangenome.pangenome_analysis import pangenome_analysis


if __name__ == "__main__":
    directory = sys.argv[1]
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('_1.fq.gz'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    tsvs = []
    target = directory + "/matrix_tsv"

    for i in [files[0]]:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        read_1, read_2 = filtering_fastq_pe(i)
        assembly = assembly_unicycler_pe(read_1, read_2)
        # assembly = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/assembly.fasta"
        bakta_annotation(assembly, name)
        annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation/" + name + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation/" + name + "_extended.tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    pangenome_analysis(target)

    asyncio.run(send_smth(text=["Ends annotation"]))
